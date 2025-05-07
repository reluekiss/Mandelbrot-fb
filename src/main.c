#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <linux/fb.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <string.h>
#include <termios.h>

#define TARGET_RE -0.743643887037158704752191506114774
#define TARGET_IM 0.131825904205311970493132056385139

#define ZOOM_SPEED 0.97f
#define BASE_ITERATIONS 100
#define ITERATION_FACTOR 100
#define THREAD_COUNT 16
#define MAX_REFERENCE_ITERATIONS 10000
#define GLITCH_THRESHOLD 1e-12

#define AUTO_ZOOM_TIMEOUT 2
#define MANUAL_ZOOM_FACTOR 0.5
#define MOVEMENT_FACTOR 0.05

#define KEY_STATE_NORMAL 0
#define KEY_STATE_ESC 1
#define KEY_STATE_BRACKET 2

double *ref_real;
double *ref_imag;
double *ref_dz_real;
double *ref_dz_imag;
int ref_iterations;

unsigned int *draw_buffer = NULL;

static const double LOG2 = 0.693147180559945;

unsigned int color_from_iteration(int iter, int max_iter, double smooth_value) {
    if (iter == max_iter) return 0;
    
    double t = (iter + 1 - smooth_value) / max_iter;
    t = t < 0 ? 0 : (t > 1 ? 1 : t); // Clamp to [0,1]
    
    unsigned char r = (unsigned char)(9 * (1-t) * t * t * t * 255);
    unsigned char g = (unsigned char)(15 * (1-t) * (1-t) * t * t * 255);
    unsigned char b = (unsigned char)(8.5 * (1-t) * (1-t) * (1-t) * t * 255);
    
    return (r << 16) | (g << 8) | b;
}

double calculate_smooth_value(double zr, double zi) {
    double mag_squared = zr*zr + zi*zi;
    return 1.0 - log(log(mag_squared) / 2.0) / LOG2;
}

void calculate_reference_orbit(double center_re, double center_im, int max_iter) {
    printf("Calculating reference orbit for (%g, %g)...\n", center_re, center_im);
    
    if (!ref_real || !ref_imag || !ref_dz_real || !ref_dz_imag) {
        // First allocation or if any buffer is NULL
        if (ref_real) free(ref_real);
        if (ref_imag) free(ref_imag);
        if (ref_dz_real) free(ref_dz_real);
        if (ref_dz_imag) free(ref_dz_imag);
        
        ref_real = malloc(max_iter * sizeof(double));
        ref_imag = malloc(max_iter * sizeof(double));
        ref_dz_real = malloc(max_iter * sizeof(double));
        ref_dz_imag = malloc(max_iter * sizeof(double));
    }
    
    double zr = 0.0, zi = 0.0;
    double dzr = 0.0, dzi = 0.0;  // Derivative values
    
    int i;
    for (i = 0; i < max_iter; i++) {
        ref_real[i] = zr;
        ref_imag[i] = zi;
        ref_dz_real[i] = dzr;
        ref_dz_imag[i] = dzi;
        
        // Check for escape
        if (zr*zr + zi*zi > 4.0)
            break;
        
        // Update derivative: dz = 2 * z * dz + 1
        double dz_temp = 2.0 * (zr * dzr - zi * dzi) + 1.0;
        dzi = 2.0 * (zr * dzi + zi * dzr);
        dzr = dz_temp;
        
        // Update z: z = z² + c
        double temp = zr*zr - zi*zi + center_re;
        zi = 2.0 * zr * zi + center_im;
        zr = temp;
    }
    
    ref_iterations = i;
    printf("Reference orbit calculated with %d iterations\n", ref_iterations);
}

typedef struct {
    struct fb_var_screeninfo *vinfo;
    struct fb_fix_screeninfo *finfo;
    double x1, y1, x2, y2;
    double center_re, center_im;
    int max_iter;
    int start_y, end_y;
    unsigned int *thread_buffer;
} ThreadData;

void* render_slice(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    int width = data->vinfo->xres;
    
    double dx = (data->x2 - data->x1) / width;
    double dy = (data->y2 - data->y1) / data->vinfo->yres;
    
    for (int y = data->start_y; y < data->end_y; y++) {
        double cy = data->y1 + y * dy;
        int buffer_y = y - data->start_y; // Relative position in thread buffer
        
        for (int x = 0; x < width; x++) {
            double cx = data->x1 + x * dx;
            
            double delta_re = cx - data->center_re;
            double delta_im = cy - data->center_im;
            
            double delta_zr = 0.0;
            double delta_zi = 0.0;
            int iter = 0;
            int escaped = 0;
            double smooth_val = 0.0;
            
            for (; iter < ref_iterations && iter < data->max_iter; iter++) {
                double ref_zr = ref_real[iter];
                double ref_zi = ref_imag[iter];
                
                // Calculate full z value (reference + delta)
                double zr = ref_zr + delta_zr;
                double zi = ref_zi + delta_zi;
                
                // Check for escape
                double mag_squared = zr*zr + zi*zi;
                if (mag_squared > 4.0) {
                    escaped = 1;
                    smooth_val = calculate_smooth_value(zr, zi);
                    break;
                }
                
                // Check for glitches (delta getting too large)
                if (fabs(delta_zr) + fabs(delta_zi) > GLITCH_THRESHOLD) {
                    zr = 0.0;
                    zi = 0.0;
                    int direct_iter;
                    for (direct_iter = 0; direct_iter < data->max_iter; direct_iter++) {
                        if (zr*zr + zi*zi > 4.0) {
                            escaped = 1;
                            smooth_val = calculate_smooth_value(zr, zi);
                            iter = direct_iter;
                            break;
                        }
                        double temp = zr*zr - zi*zi + cx;
                        zi = 2.0 * zr * zi + cy;
                        zr = temp;
                    }
                    if (!escaped) {
                        iter = data->max_iter;
                    }
                    break;
                }
                
                // Update delta using the series approximation formula
                // δz_{n+1} = 2*z_n*δz_n + δz_n² + δc
                double temp = 2.0 * (ref_zr * delta_zr - ref_zi * delta_zi) + 
                              delta_zr*delta_zr - delta_zi*delta_zi + delta_re;
                              
                delta_zi = 2.0 * (ref_zr * delta_zi + ref_zi * delta_zr) + 
                           2.0 * delta_zr * delta_zi + delta_im;
                           
                delta_zr = temp;
            }
            
            unsigned int pixel_color;
            if (!escaped && iter == data->max_iter) {
                pixel_color = 0;  // Interior point
            } else {
                pixel_color = color_from_iteration(iter, data->max_iter, smooth_val);
            }
            
            data->thread_buffer[buffer_y * width + x] = pixel_color;
        }
    }
    
    return NULL;
}

void combine_thread_buffers(ThreadData* thread_data, int thread_count) {
    struct fb_var_screeninfo *vinfo = thread_data[0].vinfo;
    int width = vinfo->xres;
    
    for (int i = 0; i < thread_count; i++) {
        ThreadData* data = &thread_data[i];
        int buffer_height = data->end_y - data->start_y;
        
        for (int y = 0; y < buffer_height; y++) {
            int src_offset = y * width;
            int dst_offset = (y + data->start_y) * width;
            
            memcpy(&draw_buffer[dst_offset], 
                   &data->thread_buffer[src_offset], 
                   width * sizeof(unsigned int));
        }
    }
}

void update_framebuffer(char *fbp, struct fb_var_screeninfo *vinfo, struct fb_fix_screeninfo *finfo) {
    int width = vinfo->xres;
    int height = vinfo->yres;
    int bytes_per_pixel = vinfo->bits_per_pixel / 8;
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            unsigned int pixel_color = draw_buffer[y * width + x];
            
            long location = (x + vinfo->xoffset) * bytes_per_pixel +
                           (y + vinfo->yoffset) * finfo->line_length;
            
            if (vinfo->bits_per_pixel == 32) {
                *(fbp + location) = pixel_color & 0xFF;             // Blue
                *(fbp + location + 1) = (pixel_color >> 8) & 0xFF;  // Green
                *(fbp + location + 2) = (pixel_color >> 16) & 0xFF; // Red
                *(fbp + location + 3) = 0;                          // Alpha
            } else if (vinfo->bits_per_pixel == 16) {
                unsigned short pixel16 = ((pixel_color >> 8) & 0xF800) |
                                         ((pixel_color >> 5) & 0x07E0) |
                                         ((pixel_color >> 3) & 0x001F);
                *((unsigned short*)(fbp + location)) = pixel16;
            } else {
                *(fbp + location) = (pixel_color & 0xFF);
            }
        }
    }
}

int process_input(double *targetX, double *targetY, double *width, double *height, 
                 int *recalculate_reference, int *auto_zoom, time_t *last_input_time) {
    static int key_state = KEY_STATE_NORMAL;
    char c;
    int moved = 0;
    
    double move_x = (*width) * MOVEMENT_FACTOR;
    double move_y = (*height) * MOVEMENT_FACTOR;
    
    while (read(STDIN_FILENO, &c, 1) > 0) {
        *last_input_time = time(NULL);
        
        if (key_state == KEY_STATE_NORMAL) {
            if (c == '\n') {
                return 1; // Exit program
            } else if (c == 27) { // ESC
                key_state = KEY_STATE_ESC;
            } else if (c == 'z' || c == 'Z') {
                *width *= MANUAL_ZOOM_FACTOR;
                *height *= MANUAL_ZOOM_FACTOR;
                moved = 1;
            } else if (c == 'x' || c == 'X') {
                *width /= MANUAL_ZOOM_FACTOR;
                *height /= MANUAL_ZOOM_FACTOR;
                moved = 1;
            } else if (c == 'c' || c == 'C') {
                *auto_zoom = *auto_zoom == 1 ? 0 : 1;
            }
        } else if (key_state == KEY_STATE_ESC) {
            if (c == '[') {
                key_state = KEY_STATE_BRACKET;
            } else {
                key_state = KEY_STATE_NORMAL; // Invalid escape sequence
            }
        } else if (key_state == KEY_STATE_BRACKET) {
            switch (c) {
                case 'A': // Up arrow
                    *targetY -= move_y;
                    moved = 1;
                    break;
                case 'B': // Down arrow
                    *targetY += move_y;
                    moved = 1;
                    break;
                case 'C': // Right arrow
                    *targetX += move_x;
                    moved = 1;
                    break;
                case 'D': // Left arrow
                    *targetX -= move_x;
                    moved = 1;
                    break;
            }
            key_state = KEY_STATE_NORMAL;
        }
    }
    
    if (moved) {
        *recalculate_reference = 1;
    }
    
    return 0;
}

int main() {
    struct termios old_tio, new_tio;
    tcgetattr(STDIN_FILENO, &old_tio);
    new_tio = old_tio;
    new_tio.c_lflag &= ~(ICANON | ECHO); // Disable canonical mode and echo
    tcsetattr(STDIN_FILENO, TCSANOW, &new_tio);
    
    int fb_fd = open("/dev/fb0", O_RDWR);
    if (fb_fd == -1) {
        perror("Error opening framebuffer device");
        tcsetattr(STDIN_FILENO, TCSANOW, &old_tio); // Restore terminal
        return 1;
    }
    
    struct fb_var_screeninfo vinfo;
    struct fb_fix_screeninfo finfo;
    if (ioctl(fb_fd, FBIOGET_FSCREENINFO, &finfo) == -1) {
        perror("Error reading fixed screen info");
        close(fb_fd);
        tcsetattr(STDIN_FILENO, TCSANOW, &old_tio);
        return 2;
    }
    if (ioctl(fb_fd, FBIOGET_VSCREENINFO, &vinfo) == -1) {
        perror("Error reading variable screen info");
        close(fb_fd);
        tcsetattr(STDIN_FILENO, TCSANOW, &old_tio);
        return 3;
    }
    
    printf("Screen resolution: %dx%d, %dbpp\n", vinfo.xres, vinfo.yres, vinfo.bits_per_pixel);
    
    long screen_size = vinfo.yres * finfo.line_length;
    char *fbp = (char *)mmap(0, screen_size, PROT_READ | PROT_WRITE, MAP_SHARED, fb_fd, 0);
    
    if ((long)fbp == -1) {
        perror("Error mapping framebuffer to memory");
        close(fb_fd);
        tcsetattr(STDIN_FILENO, TCSANOW, &old_tio);
        return 4;
    }
    
    draw_buffer = malloc(vinfo.xres * vinfo.yres * sizeof(unsigned int));
    if (!draw_buffer) {
        perror("Error allocating drawing buffer");
        munmap(fbp, screen_size);
        close(fb_fd);
        tcsetattr(STDIN_FILENO, TCSANOW, &old_tio);
        return 5;
    }
    
    ref_real = NULL;
    ref_imag = NULL;
    ref_dz_real = NULL;
    ref_dz_imag = NULL;
    
    calculate_reference_orbit(TARGET_RE, TARGET_IM, MAX_REFERENCE_ITERATIONS);
    
    double x1 = -2.5;
    double y1 = -1.5;
    double x2 = 1.0;
    double y2 = 1.5;
    
    double targetX = TARGET_RE;
    double targetY = TARGET_IM;
    
    pthread_t threads[THREAD_COUNT];
    ThreadData thread_data[THREAD_COUNT];
    
    int slice_height = vinfo.yres / THREAD_COUNT;
    for (int i = 0; i < THREAD_COUNT; i++) {
        int end_y = (i == THREAD_COUNT - 1) ? vinfo.yres : (i + 1) * slice_height;
        int buffer_height = end_y - (i * slice_height);
        
        thread_data[i].thread_buffer = malloc(vinfo.xres * buffer_height * sizeof(unsigned int));
        if (!thread_data[i].thread_buffer) {
            perror("Error allocating thread buffer");
            goto end;
        }
    }
    
    int flags = fcntl(STDIN_FILENO, F_GETFL, 0);
    fcntl(STDIN_FILENO, F_SETFL, flags | O_NONBLOCK);

int auto_zoom = 1; // Start with auto-zoom enabled
    time_t last_input_time = time(NULL);
    
    printf("Mandelbrot Explorer with Navigation\n");
    printf("Use arrow keys to move around, z to zoom in, x to zoom out, c to toggle auto zoom, press Enter to exit.\n");
    printf("Auto-zoom will resume after 2 seconds of inactivity.\n");
    
    int running = 1;
    int frame_count = 0;
    time_t start_time = time(NULL);
    double initial_width = 3.5;
    int recalculate_reference = 0;
    
    while (running) {
        double width = x2 - x1;
        double height = y2 - y1;
        
        if (process_input(&targetX, &targetY, &width, &height, &recalculate_reference, &auto_zoom, &last_input_time)) {
            running = 0;
            break;
        }
        
        x1 = targetX - width / 2;
        x2 = targetX + width / 2;
        y1 = targetY - height / 2;
        y2 = targetY + height / 2;
        
        time_t current_time = time(NULL);
        if (!auto_zoom && difftime(current_time, last_input_time) >= AUTO_ZOOM_TIMEOUT) {
            auto_zoom = 1;
            printf("\rAuto-zoom resumed                                      \n");
        }
        
        if (recalculate_reference) {
            calculate_reference_orbit(targetX, targetY, MAX_REFERENCE_ITERATIONS);
            recalculate_reference = 0;
        }
        
        double zoom_level = initial_width / width;
        
        int max_iter = BASE_ITERATIONS + (int)(ITERATION_FACTOR * log10(zoom_level));
        if (max_iter < BASE_ITERATIONS) max_iter = BASE_ITERATIONS;
        if (max_iter > MAX_REFERENCE_ITERATIONS) max_iter = MAX_REFERENCE_ITERATIONS;
        
        for (int i = 0; i < THREAD_COUNT; i++) {
            thread_data[i].vinfo = &vinfo;
            thread_data[i].finfo = &finfo;
            thread_data[i].x1 = x1;
            thread_data[i].y1 = y1;
            thread_data[i].x2 = x2;
            thread_data[i].y2 = y2;
            thread_data[i].center_re = targetX;
            thread_data[i].center_im = targetY;
            thread_data[i].max_iter = max_iter;
            thread_data[i].start_y = i * slice_height;
            thread_data[i].end_y = (i == THREAD_COUNT - 1) ? vinfo.yres : (i + 1) * slice_height;
            
            pthread_create(&threads[i], NULL, render_slice, &thread_data[i]);
        }
        
        for (int i = 0; i < THREAD_COUNT; i++) {
            pthread_join(threads[i], NULL);
        }
        
        combine_thread_buffers(thread_data, THREAD_COUNT);
        
        update_framebuffer(fbp, &vinfo, &finfo);
        
        if (auto_zoom && width > 1e-14) {
            double newWidth = width * ZOOM_SPEED;
            double newHeight = height * ZOOM_SPEED;
            
            x1 = targetX - newWidth / 2;
            x2 = targetX + newWidth / 2;
            y1 = targetY - newHeight / 2;
            y2 = targetY + newHeight / 2;
        } else if (width <= 1e-14) {
            // Reset if we've zoomed in too far
            x1 = -2.5;
            y1 = -1.5;
            x2 = 1.0;
            y2 = 1.5;
        }
        
        frame_count++;
        current_time = time(NULL);
        double elapsed = difftime(current_time, start_time);
        
        if (elapsed >= 1.0) {
            const char* zoom_mode = auto_zoom ? "Auto" : "Manual";
            printf("\rZoom: %.2e×   Mode: %s   Position: (%.15g, %.15g)   Iterations: %d   FPS: %.1f   ",
                   zoom_level, zoom_mode, targetX, targetY, max_iter, frame_count / elapsed);
            fflush(stdout);
            
            start_time = current_time;
            frame_count = 0;
        }
    }
    
    printf("\nExiting...\n");
    
end:
    for (int i = 0; i < THREAD_COUNT; i++) {
        if (thread_data[i].thread_buffer)
            free(thread_data[i].thread_buffer);
    }

    free(draw_buffer);
    free(ref_real);
    free(ref_imag);
    free(ref_dz_real);
    free(ref_dz_imag);
    munmap(fbp, screen_size);
    close(fb_fd);
    
    tcsetattr(STDIN_FILENO, TCSANOW, &old_tio);
    
    return 0;
}
