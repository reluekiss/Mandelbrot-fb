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

#define TARGET_RE -0.743643887037158704752191506114774
#define TARGET_IM 0.131825904205311970493132056385139

#define ZOOM_SPEED 0.97f
#define BASE_ITERATIONS 100
#define ITERATION_FACTOR 100
#define THREAD_COUNT 16

unsigned int color_from_iteration(int iter, int max_iter) {
    if (iter == max_iter) return 0;
    
    float t = (float)iter / max_iter;
    
    unsigned char r = (unsigned char)(9 * (1-t) * t * t * t * 255);
    unsigned char g = (unsigned char)(15 * (1-t) * (1-t) * t * t * 255);
    unsigned char b = (unsigned char)(8.5 * (1-t) * (1-t) * (1-t) * t * 255);
    
    return (r << 16) | (g << 8) | b;
}

typedef struct {
    char *fbp;
    struct fb_var_screeninfo *vinfo;
    struct fb_fix_screeninfo *finfo;
    double x1, y1, x2, y2;
    int max_iter;
    int start_y, end_y;
} ThreadData;

void* render_slice(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    int bytes_per_pixel = data->vinfo->bits_per_pixel / 8;
    int line_length = data->finfo->line_length;
    
    double dx = (data->x2 - data->x1) / data->vinfo->xres;
    double dy = (data->y2 - data->y1) / data->vinfo->yres;
    
    for (int y = data->start_y; y < data->end_y; y++) {
        double cy = data->y1 + y * dy;
        
        for (uint x = 0; x < data->vinfo->xres; x++) {
            double cx = data->x1 + x * dx;
            
            double zx = 0;
            double zy = 0;
            double zx2 = 0;
            double zy2 = 0;
            int iter = 0;
            
            while (iter < data->max_iter && (zx2 + zy2 < 4.0)) {
                zy = 2 * zx * zy + cy;
                zx = zx2 - zy2 + cx;
                zx2 = zx * zx;
                zy2 = zy * zy;
                iter++;
            }
            
            long location = (x + data->vinfo->xoffset) * bytes_per_pixel +
                           (y + data->vinfo->yoffset) * line_length;
            
            unsigned int pixel_color = color_from_iteration(iter, data->max_iter);
            
            if (data->vinfo->bits_per_pixel == 32) {
                *(data->fbp + location) = pixel_color & 0xFF;               // Blue
                *(data->fbp + location + 1) = (pixel_color >> 8) & 0xFF;    // Green
                *(data->fbp + location + 2) = (pixel_color >> 16) & 0xFF;   // Red
                *(data->fbp + location + 3) = 0;                            // Alpha (ignored)
            } else if (data->vinfo->bits_per_pixel == 16) {
                unsigned short pixel16 = ((pixel_color >> 8) & 0xF800) |
                                         ((pixel_color >> 5) & 0x07E0) |
                                         ((pixel_color >> 3) & 0x001F);
                *((unsigned short*)(data->fbp + location)) = pixel16;
            } else {
                *(data->fbp + location) = (iter * 255) / data->max_iter;
            }
        }
    }
    
    return NULL;
}

int main() {
    int fb_fd = open("/dev/fb0", O_RDWR);
    if (fb_fd == -1) {
        perror("Error opening framebuffer device");
        return 1;
    }
    
    struct fb_var_screeninfo vinfo;
    struct fb_fix_screeninfo finfo;
    if (ioctl(fb_fd, FBIOGET_FSCREENINFO, &finfo) == -1) {
        perror("Error reading fixed screen info");
        close(fb_fd);
        return 2;
    }
    if (ioctl(fb_fd, FBIOGET_VSCREENINFO, &vinfo) == -1) {
        perror("Error reading variable screen info");
        close(fb_fd);
        return 3;
    }
    
    printf("Screen resolution: %dx%d, %dbpp\n", vinfo.xres, vinfo.yres, vinfo.bits_per_pixel);
    
    long screen_size = vinfo.yres * finfo.line_length;
    char *fbp = (char *)mmap(0, screen_size, PROT_READ | PROT_WRITE, MAP_SHARED, fb_fd, 0);
    
    if ((long)fbp == -1) {
        perror("Error mapping framebuffer to memory");
        close(fb_fd);
        return 4;
    }
    
    double x1 = -2.5;
    double y1 = -1.5;
    double x2 = 1.0;
    double y2 = 1.5;
    
    double targetX = TARGET_RE;
    double targetY = TARGET_IM;
    
    pthread_t threads[THREAD_COUNT];
    ThreadData thread_data[THREAD_COUNT];
    
    int flags = fcntl(STDIN_FILENO, F_GETFL, 0);
    fcntl(STDIN_FILENO, F_SETFL, flags | O_NONBLOCK);
    
    printf("Rendering Mandelbrot set...\n");
    printf("Press Enter to exit.\n");
    
    int running = 1;
    int frame_count = 0;
    time_t start_time = time(NULL);
    double initial_width = 3.5;
    
    while (running) {
        char c;
        if (read(STDIN_FILENO, &c, 1) > 0) {
            if (c == '\n') {
                running = 0;
                break;
            }
        }
        
        double width = x2 - x1;
        double height = y2 - y1;
        
        double zoom_level = initial_width / width;
        
        int max_iter = BASE_ITERATIONS + (int)(ITERATION_FACTOR * log10(zoom_level));
        if (max_iter < BASE_ITERATIONS) max_iter = BASE_ITERATIONS;
        
        int slice_height = vinfo.yres / THREAD_COUNT;
        
        for (int i = 0; i < THREAD_COUNT; i++) {
            thread_data[i].fbp = fbp;
            thread_data[i].vinfo = &vinfo;
            thread_data[i].finfo = &finfo;
            thread_data[i].x1 = x1;
            thread_data[i].y1 = y1;
            thread_data[i].x2 = x2;
            thread_data[i].y2 = y2;
            thread_data[i].max_iter = max_iter;
            thread_data[i].start_y = i * slice_height;
            thread_data[i].end_y = (i == THREAD_COUNT - 1) ? vinfo.yres : (i + 1) * slice_height;
            
            pthread_create(&threads[i], NULL, render_slice, &thread_data[i]);
        }
        
        for (int i = 0; i < THREAD_COUNT; i++) {
            pthread_join(threads[i], NULL);
        }
        
        frame_count++;
        time_t current_time = time(NULL);
        double elapsed = difftime(current_time, start_time);
        
        if (elapsed >= 1.0) {
            printf("\rZoom: %.0f×   Iterations: %d   FPS: %.1f   ",
                   zoom_level, max_iter, frame_count / elapsed);
            fflush(stdout);
            
            start_time = current_time;
            frame_count = 0;
        }
        
        if (width > 1e-14) {
            double newWidth = width * ZOOM_SPEED;
            double newHeight = height * ZOOM_SPEED;
            
            x1 = targetX - newWidth / 2;
            x2 = targetX + newWidth / 2;
            y1 = targetY - newHeight / 2;
            y2 = targetY + newHeight / 2;
        } else {
            x1 = -2.5;
            y1 = -1.5;
            x2 = 1.0;
            y2 = 1.5;
        }
    }
    
    printf("\nExiting...\n");
    
    munmap(fbp, screen_size);
    close(fb_fd);
    
    return 0;
}
