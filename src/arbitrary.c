#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <pthread.h>
#include <math.h>
#include <linux/fb.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <string.h>
#include <termios.h>
#include <mpfr.h>

#define PRECISION 256

#define TARGET_RE "-0.743643887037158704752191506114774"
#define TARGET_IM "0.131825904205311970493132056385139"

#define ZOOM_SPEED "0.97"
#define BASE_ITERATIONS 100
#define ITERATION_FACTOR 100
#define THREAD_COUNT 16
#define MAX_REFERENCE_ITERATIONS 10000
#define GLITCH_THRESHOLD "1e-12"

#define AUTO_ZOOM_TIMEOUT 2
#define MANUAL_ZOOM_FACTOR "0.5"
#define MOVEMENT_FACTOR "0.05"

#define KEY_STATE_NORMAL 0
#define KEY_STATE_ESC 1
#define KEY_STATE_BRACKET 2

mpfr_t *ref_real;
mpfr_t *ref_imag;
mpfr_t *ref_dz_real;
mpfr_t *ref_dz_imag;
int ref_iterations;

unsigned int *draw_buffer = NULL;

static const double LOG2 = 0.693147180559945;

static inline unsigned int color_from_iteration(int iter, int max_iter, double smooth_value) {
    if (iter == max_iter) return 0;
    double t = (iter + 1 - smooth_value) / max_iter;
    t = t < 0 ? 0 : (t > 1 ? 1 : t); // Clamp to [0,1]
    unsigned char r = (unsigned char)(9 * (1-t) * t * t * t * 255);
    unsigned char g = (unsigned char)(15 * (1-t) * (1-t) * t * t * 255);
    unsigned char b = (unsigned char)(8.5 * (1-t) * (1-t) * (1-t) * t * 255);
    return (r << 16) | (g << 8) | b;
}

static inline double calculate_smooth_value(mpfr_t zr, mpfr_t zi) {
    mpfr_t mag_squared, log_mag, log_log_mag, result;
    mpfr_inits2(PRECISION, mag_squared, log_mag, log_log_mag, result, (mpfr_ptr)0);

    mpfr_sqr(mag_squared, zr, MPFR_RNDN);
    mpfr_t temp;
    mpfr_init2(temp, PRECISION);
    mpfr_sqr(temp, zi, MPFR_RNDN);
    mpfr_add(mag_squared, mag_squared, temp, MPFR_RNDN);

    mpfr_log(log_mag, mag_squared, MPFR_RNDN);
    mpfr_div_d(log_mag, log_mag, 2.0, MPFR_RNDN);
    mpfr_log(log_log_mag, log_mag, MPFR_RNDN);
    mpfr_set_d(result, LOG2, MPFR_RNDN);
    mpfr_div(log_log_mag, log_log_mag, result, MPFR_RNDN);
    mpfr_ui_sub(result, 1.0, log_log_mag, MPFR_RNDN);

    double ret = mpfr_get_d(result, MPFR_RNDN);

    mpfr_clears(mag_squared, log_mag, log_log_mag, result, temp, (mpfr_ptr)0);
    return ret;
}

static inline void calculate_reference_orbit(mpfr_t center_re, mpfr_t center_im, int max_iter) {
    printf("Calculating reference orbit...\n");

    if (ref_real) {
        for (int i = 0; i < MAX_REFERENCE_ITERATIONS; i++) {
            mpfr_clear(ref_real[i]);
            mpfr_clear(ref_imag[i]);
            mpfr_clear(ref_dz_real[i]);
            mpfr_clear(ref_dz_imag[i]);
        }
        free(ref_real); free(ref_imag); free(ref_dz_real); free(ref_dz_imag);
    }

    ref_real = malloc(max_iter * sizeof(mpfr_t));
    ref_imag = malloc(max_iter * sizeof(mpfr_t));
    ref_dz_real = malloc(max_iter * sizeof(mpfr_t));
    ref_dz_imag = malloc(max_iter * sizeof(mpfr_t));
    for (int i = 0; i < max_iter; i++) {
        mpfr_init2(ref_real[i], PRECISION);
        mpfr_init2(ref_imag[i], PRECISION);
        mpfr_init2(ref_dz_real[i], PRECISION);
        mpfr_init2(ref_dz_imag[i], PRECISION);
    }

    mpfr_t zr, zi, dzr, dzi, temp, dz_temp, mag_squared;
    mpfr_inits2(PRECISION, zr, zi, dzr, dzi, temp, dz_temp, mag_squared, (mpfr_ptr)0);

    mpfr_set_d(zr, 0.0, MPFR_RNDN);
    mpfr_set_d(zi, 0.0, MPFR_RNDN);
    mpfr_set_d(dzr, 0.0, MPFR_RNDN);
    mpfr_set_d(dzi, 0.0, MPFR_RNDN);

    int i;
    for (i = 0; i < max_iter; i++) {
        mpfr_set(ref_real[i], zr, MPFR_RNDN);
        mpfr_set(ref_imag[i], zi, MPFR_RNDN);
        mpfr_set(ref_dz_real[i], dzr, MPFR_RNDN);
        mpfr_set(ref_dz_imag[i], dzi, MPFR_RNDN);

        // mag_squared = zr*zr + zi*zi
        mpfr_sqr(mag_squared, zr, MPFR_RNDN);
        mpfr_sqr(temp, zi, MPFR_RNDN);
        mpfr_add(mag_squared, mag_squared, temp, MPFR_RNDN);

        if (mpfr_cmp_d(mag_squared, 4.0) > 0)
            break;

        // dz_temp = 2.0 * (zr * dzr - zi * dzi) + 1.0;
        mpfr_mul(temp, zr, dzr, MPFR_RNDN);
        mpfr_mul(dz_temp, zi, dzi, MPFR_RNDN);
        mpfr_sub(temp, temp, dz_temp, MPFR_RNDN);
        mpfr_mul_d(temp, temp, 2.0, MPFR_RNDN);
        mpfr_add_d(dz_temp, temp, 1.0, MPFR_RNDN);

        // dzi = 2.0 * (zr * dzi + zi * dzr);
        mpfr_mul(temp, zr, dzi, MPFR_RNDN);
        mpfr_mul(dzi, zi, dzr, MPFR_RNDN);
        mpfr_add(temp, temp, dzi, MPFR_RNDN);
        mpfr_mul_d(dzi, temp, 2.0, MPFR_RNDN);

        // dzr = dz_temp;
        mpfr_set(dzr, dz_temp, MPFR_RNDN);

        // temp = zr*zr - zi*zi + center_re;
        mpfr_sqr(temp, zr, MPFR_RNDN);
        mpfr_sqr(dz_temp, zi, MPFR_RNDN);
        mpfr_sub(temp, temp, dz_temp, MPFR_RNDN);
        mpfr_add(temp, temp, center_re, MPFR_RNDN);

        // zi = 2.0 * zr * zi + center_im;
        mpfr_mul(dz_temp, zr, zi, MPFR_RNDN);
        mpfr_mul_d(dz_temp, dz_temp, 2.0, MPFR_RNDN);
        mpfr_add(zi, dz_temp, center_im, MPFR_RNDN);

        // zr = temp;
        mpfr_set(zr, temp, MPFR_RNDN);
    }

    ref_iterations = i;

    mpfr_clears(zr, zi, dzr, dzi, temp, dz_temp, mag_squared, (mpfr_ptr)0);
}

typedef struct {
    struct fb_var_screeninfo *vinfo;
    struct fb_fix_screeninfo *finfo;
    mpfr_t x1, y1, x2, y2;
    mpfr_t center_re, center_im;
    int max_iter;
    int start_y, end_y;
    unsigned int *thread_buffer;
} ThreadData;

static inline void* render_slice(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    int width = data->vinfo->xres;
    int height = data->vinfo->yres;

    mpfr_t dx, dy, cx, cy, delta_re, delta_im, delta_zr, delta_zi;
    mpfr_t ref_zr, ref_zi, zr, zi, mag_squared, temp, temp2, glitch_threshold;
    mpfr_inits2(PRECISION, dx, dy, cx, cy, delta_re, delta_im, delta_zr, delta_zi,
                ref_zr, ref_zi, zr, zi, mag_squared, temp, temp2, glitch_threshold, (mpfr_ptr)0);

    mpfr_sub(temp, data->x2, data->x1, MPFR_RNDN);
    mpfr_set_si(temp2, width, MPFR_RNDN);
    mpfr_div(dx, temp, temp2, MPFR_RNDN);

    mpfr_sub(temp, data->y2, data->y1, MPFR_RNDN);
    mpfr_set_si(temp2, height, MPFR_RNDN);
    mpfr_div(dy, temp, temp2, MPFR_RNDN);

    mpfr_set_str(glitch_threshold, GLITCH_THRESHOLD, 10, MPFR_RNDN);

    for (int y = data->start_y; y < data->end_y; y++) {
        mpfr_set_si(temp, y, MPFR_RNDN);
        mpfr_mul(cy, temp, dy, MPFR_RNDN);
        mpfr_add(cy, data->y1, cy, MPFR_RNDN);

        int buffer_y = y - data->start_y;
        for (int x = 0; x < width; x++) {
            mpfr_set_si(temp, x, MPFR_RNDN);
            mpfr_mul(cx, temp, dx, MPFR_RNDN);
            mpfr_add(cx, data->x1, cx, MPFR_RNDN);

            mpfr_sub(delta_re, cx, data->center_re, MPFR_RNDN);
            mpfr_sub(delta_im, cy, data->center_im, MPFR_RNDN);

            mpfr_set_d(delta_zr, 0.0, MPFR_RNDN);
            mpfr_set_d(delta_zi, 0.0, MPFR_RNDN);

            int iter = 0;
            int escaped = 0;
            double smooth_val = 0.0;

            for (; iter < ref_iterations && iter < data->max_iter; iter++) {
                mpfr_set(ref_zr, ref_real[iter], MPFR_RNDN);
                mpfr_set(ref_zi, ref_imag[iter], MPFR_RNDN);

                mpfr_add(zr, ref_zr, delta_zr, MPFR_RNDN);
                mpfr_add(zi, ref_zi, delta_zi, MPFR_RNDN);

                mpfr_sqr(mag_squared, zr, MPFR_RNDN);
                mpfr_sqr(temp, zi, MPFR_RNDN);
                mpfr_add(mag_squared, mag_squared, temp, MPFR_RNDN);

                if (mpfr_cmp_d(mag_squared, 4.0) > 0) {
                    escaped = 1;
                    smooth_val = calculate_smooth_value(zr, zi);
                    break;
                }

                mpfr_abs(temp, delta_zr, MPFR_RNDN);
                mpfr_abs(temp2, delta_zi, MPFR_RNDN);
                mpfr_add(temp, temp, temp2, MPFR_RNDN);
                if (mpfr_cmp(temp, glitch_threshold) > 0) {
                    // Fallback to direct calculation
                    mpfr_set_d(zr, 0.0, MPFR_RNDN);
                    mpfr_set_d(zi, 0.0, MPFR_RNDN);
                    int direct_iter;
                    for (direct_iter = 0; direct_iter < data->max_iter; direct_iter++) {
                        mpfr_sqr(mag_squared, zr, MPFR_RNDN);
                        mpfr_sqr(temp, zi, MPFR_RNDN);
                        mpfr_add(mag_squared, mag_squared, temp, MPFR_RNDN);
                        if (mpfr_cmp_d(mag_squared, 4.0) > 0) {
                            escaped = 1;
                            smooth_val = calculate_smooth_value(zr, zi);
                            iter = direct_iter;
                            break;
                        }
                        // temp = zr*zr - zi*zi + cx;
                        mpfr_sqr(temp, zr, MPFR_RNDN);
                        mpfr_sqr(temp2, zi, MPFR_RNDN);
                        mpfr_sub(temp, temp, temp2, MPFR_RNDN);
                        mpfr_add(temp, temp, cx, MPFR_RNDN);

                        // zi = 2.0 * zr * zi + cy;
                        mpfr_mul(temp2, zr, zi, MPFR_RNDN);
                        mpfr_mul_d(temp2, temp2, 2.0, MPFR_RNDN);
                        mpfr_add(zi, temp2, cy, MPFR_RNDN);

                        mpfr_set(zr, temp, MPFR_RNDN);
                    }
                    if (!escaped) iter = data->max_iter;
                    break;
                }

                // delta update: δz_{n+1} = 2*z_n*δz_n + δz_n² + δc
                // temp = 2.0 * (ref_zr * delta_zr - ref_zi * delta_zi)
                mpfr_mul(temp, ref_zr, delta_zr, MPFR_RNDN);
                mpfr_mul(temp2, ref_zi, delta_zi, MPFR_RNDN);
                mpfr_sub(temp, temp, temp2, MPFR_RNDN);
                mpfr_mul_d(temp, temp, 2.0, MPFR_RNDN);

                // + delta_zr*delta_zr - delta_zi*delta_zi
                mpfr_sqr(temp2, delta_zr, MPFR_RNDN);
                mpfr_add(temp, temp, temp2, MPFR_RNDN);
                mpfr_sqr(temp2, delta_zi, MPFR_RNDN);
                mpfr_sub(temp, temp, temp2, MPFR_RNDN);

                // + delta_re
                mpfr_add(temp, temp, delta_re, MPFR_RNDN);

                // delta_zi = 2.0 * (ref_zr * delta_zi + ref_zi * delta_zr) + 2.0 * delta_zr * delta_zi + delta_im;
                mpfr_mul(temp2, ref_zr, delta_zi, MPFR_RNDN);
                mpfr_mul(mag_squared, ref_zi, delta_zr, MPFR_RNDN);
                mpfr_add(temp2, temp2, mag_squared, MPFR_RNDN);
                mpfr_mul_d(temp2, temp2, 2.0, MPFR_RNDN);

                mpfr_mul(mag_squared, delta_zr, delta_zi, MPFR_RNDN);
                mpfr_mul_d(mag_squared, mag_squared, 2.0, MPFR_RNDN);
                mpfr_add(temp2, temp2, mag_squared, MPFR_RNDN);
                mpfr_add(temp2, temp2, delta_im, MPFR_RNDN);

                mpfr_set(delta_zi, temp2, MPFR_RNDN);
                mpfr_set(delta_zr, temp, MPFR_RNDN);
            }

            unsigned int pixel_color;
            if (!escaped && iter == data->max_iter) {
                pixel_color = 0;
            } else {
                pixel_color = color_from_iteration(iter, data->max_iter, smooth_val);
            }
            data->thread_buffer[buffer_y * width + x] = pixel_color;
        }
    }

    mpfr_clears(dx, dy, cx, cy, delta_re, delta_im, delta_zr, delta_zi,
                ref_zr, ref_zi, zr, zi, mag_squared, temp, temp2, glitch_threshold, (mpfr_ptr)0);
    return NULL;
}

static inline void combine_thread_buffers(ThreadData* thread_data, int thread_count) {
    struct fb_var_screeninfo *vinfo = thread_data[0].vinfo;
    int width = vinfo->xres;
    for (int i = 0; i < thread_count; i++) {
        ThreadData* data = &thread_data[i];
        int buffer_height = data->end_y - data->start_y;
        for (int y = 0; y < buffer_height; y++) {
            int src_offset = y * width;
            int dst_offset = (y + data->start_y) * width;
            memcpy(&draw_buffer[dst_offset], &data->thread_buffer[src_offset], width * sizeof(unsigned int));
        }
    }
}

static inline void update_framebuffer(char *fbp, struct fb_var_screeninfo *vinfo, struct fb_fix_screeninfo *finfo) {
    int width = vinfo->xres;
    int height = vinfo->yres;
    int bytes_per_pixel = vinfo->bits_per_pixel / 8;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            unsigned int pixel_color = draw_buffer[y * width + x];
            long location = (x + vinfo->xoffset) * bytes_per_pixel +
                            (y + vinfo->yoffset) * finfo->line_length;
            if (vinfo->bits_per_pixel == 32) {
                *(fbp + location) = pixel_color & 0xFF;
                *(fbp + location + 1) = (pixel_color >> 8) & 0xFF;
                *(fbp + location + 2) = (pixel_color >> 16) & 0xFF;
                *(fbp + location + 3) = 0;
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

static inline int process_input(mpfr_t targetX, mpfr_t targetY, mpfr_t width, mpfr_t height,
                               int *recalculate_reference, int *auto_zoom, time_t *last_input_time) {
    static int key_state = KEY_STATE_NORMAL;
    char c;
    int moved = 0;

    mpfr_t move_x, move_y, temp;
    mpfr_inits2(PRECISION, move_x, move_y, temp, (mpfr_ptr)0);

    mpfr_set_str(temp, MOVEMENT_FACTOR, 10, MPFR_RNDN);
    mpfr_mul(move_x, width, temp, MPFR_RNDN);
    mpfr_mul(move_y, height, temp, MPFR_RNDN);

    while (read(STDIN_FILENO, &c, 1) > 0) {
        *last_input_time = time(NULL);
        if (key_state == KEY_STATE_NORMAL) {
            if (c == '\n') {
                mpfr_clears(move_x, move_y, temp, (mpfr_ptr)0);
                return 1;
            } else if (c == 27) {
                key_state = KEY_STATE_ESC;
            } else if (c == 'z' || c == 'Z') {
                mpfr_set_str(temp, MANUAL_ZOOM_FACTOR, 10, MPFR_RNDN);
                mpfr_mul(width, width, temp, MPFR_RNDN);
                mpfr_mul(height, height, temp, MPFR_RNDN);
                moved = 1;
            } else if (c == 'x' || c == 'X') {
                mpfr_set_str(temp, MANUAL_ZOOM_FACTOR, 10, MPFR_RNDN);
                mpfr_div(width, width, temp, MPFR_RNDN);
                mpfr_div(height, height, temp, MPFR_RNDN);
                moved = 1;
            } else if (c == 'c' || c == 'C') {
                *auto_zoom = *auto_zoom == 1 ? 0 : 1;
            }
        } else if (key_state == KEY_STATE_ESC) {
            if (c == '[') {
                key_state = KEY_STATE_BRACKET;
            } else {
                key_state = KEY_STATE_NORMAL;
            }
        } else if (key_state == KEY_STATE_BRACKET) {
            switch (c) {
                case 'A':
                    mpfr_sub(targetY, targetY, move_y, MPFR_RNDN);
                    moved = 1;
                    break;
                case 'B':
                    mpfr_add(targetY, targetY, move_y, MPFR_RNDN);
                    moved = 1;
                    break;
                case 'C':
                    mpfr_add(targetX, targetX, move_x, MPFR_RNDN);
                    moved = 1;
                    break;
                case 'D':
                    mpfr_sub(targetX, targetX, move_x, MPFR_RNDN);
                    moved = 1;
                    break;
            }
            key_state = KEY_STATE_NORMAL;
        }
    }
    if (moved) *recalculate_reference = 1;
    mpfr_clears(move_x, move_y, temp, (mpfr_ptr)0);
    return 0;
}

int main() {
    struct termios old_tio, new_tio;
    tcgetattr(STDIN_FILENO, &old_tio);
    new_tio = old_tio;
    new_tio.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &new_tio);

    int fb_fd = open("/dev/fb0", O_RDWR);
    if (fb_fd == -1) {
        perror("Error opening framebuffer device");
        tcsetattr(STDIN_FILENO, TCSANOW, &old_tio);
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

    mpfr_t targetX, targetY, x1, y1, x2, y2, width, height, initial_width, temp;
    mpfr_inits2(PRECISION, targetX, targetY, x1, y1, x2, y2, width, height, initial_width, temp, (mpfr_ptr)0);

    mpfr_set_str(targetX, TARGET_RE, 10, MPFR_RNDN);
    mpfr_set_str(targetY, TARGET_IM, 10, MPFR_RNDN);

    mpfr_set_d(x1, -2.5, MPFR_RNDN);
    mpfr_set_d(y1, -1.5, MPFR_RNDN);
    mpfr_set_d(x2, 1.0, MPFR_RNDN);
    mpfr_set_d(y2, 1.5, MPFR_RNDN);

    mpfr_sub(temp, x2, x1, MPFR_RNDN);
    mpfr_set(initial_width, temp, MPFR_RNDN);

    mpfr_t ref_center_re, ref_center_im;
    mpfr_inits2(PRECISION, ref_center_re, ref_center_im, (mpfr_ptr)0);
    mpfr_set(ref_center_re, targetX, MPFR_RNDN);
    mpfr_set(ref_center_im, targetY, MPFR_RNDN);

    calculate_reference_orbit(ref_center_re, ref_center_im, MAX_REFERENCE_ITERATIONS);

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

    int auto_zoom = 1;
    time_t last_input_time = time(NULL);

    printf("Mandelbrot Explorer (Arbitrary Precision)\n");
    printf("Use arrow keys to move, z to zoom in, x to zoom out, c to toggle auto zoom, Enter to exit.\n");

    int running = 1;
    int frame_count = 0;
    time_t start_time = time(NULL);
    int recalculate_reference = 0;

    mpfr_t zoom_level, zoom_speed, min_width, manual_zoom_factor;
    mpfr_inits2(PRECISION, zoom_level, zoom_speed, min_width, manual_zoom_factor, (mpfr_ptr)0);
    mpfr_set_str(zoom_speed, ZOOM_SPEED, 10, MPFR_RNDN);
    mpfr_set_str(min_width, "1e-14", 10, MPFR_RNDN);
    mpfr_set_str(manual_zoom_factor, MANUAL_ZOOM_FACTOR, 10, MPFR_RNDN);

    while (running) {
        mpfr_sub(width, x2, x1, MPFR_RNDN);
        mpfr_sub(height, y2, y1, MPFR_RNDN);

        if (process_input(targetX, targetY, width, height, &recalculate_reference, &auto_zoom, &last_input_time)) {
            running = 0;
            break;
        }

        mpfr_div_ui(temp, width, 2, MPFR_RNDN);
        mpfr_sub(x1, targetX, temp, MPFR_RNDN);
        mpfr_add(x2, targetX, temp, MPFR_RNDN);

        mpfr_div_ui(temp, height, 2, MPFR_RNDN);
        mpfr_sub(y1, targetY, temp, MPFR_RNDN);
        mpfr_add(y2, targetY, temp, MPFR_RNDN);

        time_t current_time = time(NULL);
        if (!auto_zoom && difftime(current_time, last_input_time) >= AUTO_ZOOM_TIMEOUT) {
            auto_zoom = 1;
            printf("\rAuto-zoom resumed                                      \n");
        }

        if (recalculate_reference) {
            mpfr_set(ref_center_re, targetX, MPFR_RNDN);
            mpfr_set(ref_center_im, targetY, MPFR_RNDN);
            calculate_reference_orbit(ref_center_re, ref_center_im, MAX_REFERENCE_ITERATIONS);
            recalculate_reference = 0;
        }

        mpfr_div(zoom_level, initial_width, width, MPFR_RNDN);
        double zoom_level_d = mpfr_get_d(zoom_level, MPFR_RNDN);

        int max_iter = BASE_ITERATIONS + (int)(ITERATION_FACTOR * log10(zoom_level_d));
        if (max_iter < BASE_ITERATIONS) max_iter = BASE_ITERATIONS;
        if (max_iter > MAX_REFERENCE_ITERATIONS) max_iter = MAX_REFERENCE_ITERATIONS;

        for (int i = 0; i < THREAD_COUNT; i++) {
            thread_data[i].vinfo = &vinfo;
            thread_data[i].finfo = &finfo;
            mpfr_init2(thread_data[i].x1, PRECISION);
            mpfr_init2(thread_data[i].y1, PRECISION);
            mpfr_init2(thread_data[i].x2, PRECISION);
            mpfr_init2(thread_data[i].y2, PRECISION);
            mpfr_init2(thread_data[i].center_re, PRECISION);
            mpfr_init2(thread_data[i].center_im, PRECISION);

            mpfr_set(thread_data[i].x1, x1, MPFR_RNDN);
            mpfr_set(thread_data[i].y1, y1, MPFR_RNDN);
            mpfr_set(thread_data[i].x2, x2, MPFR_RNDN);
            mpfr_set(thread_data[i].y2, y2, MPFR_RNDN);
            mpfr_set(thread_data[i].center_re, targetX, MPFR_RNDN);
            mpfr_set(thread_data[i].center_im, targetY, MPFR_RNDN);

            thread_data[i].max_iter = max_iter;
            thread_data[i].start_y = i * slice_height;
            thread_data[i].end_y = (i == THREAD_COUNT - 1) ? vinfo.yres : (i + 1) * slice_height;

            pthread_create(&threads[i], NULL, render_slice, &thread_data[i]);
        }

        for (int i = 0; i < THREAD_COUNT; i++) {
            pthread_join(threads[i], NULL);
            mpfr_clear(thread_data[i].x1);
            mpfr_clear(thread_data[i].y1);
            mpfr_clear(thread_data[i].x2);
            mpfr_clear(thread_data[i].y2);
            mpfr_clear(thread_data[i].center_re);
            mpfr_clear(thread_data[i].center_im);
        }

        combine_thread_buffers(thread_data, THREAD_COUNT);
        update_framebuffer(fbp, &vinfo, &finfo);

        if (auto_zoom && mpfr_cmp(width, min_width) > 0) {
            mpfr_mul(width, width, zoom_speed, MPFR_RNDN);
            mpfr_mul(height, height, zoom_speed, MPFR_RNDN);

            mpfr_div_ui(temp, width, 2, MPFR_RNDN);
            mpfr_sub(x1, targetX, temp, MPFR_RNDN);
            mpfr_add(x2, targetX, temp, MPFR_RNDN);

            mpfr_div_ui(temp, height, 2, MPFR_RNDN);
            mpfr_sub(y1, targetY, temp, MPFR_RNDN);
            mpfr_add(y2, targetY, temp, MPFR_RNDN);
        } else if (mpfr_cmp(width, min_width) <= 0) {
            mpfr_set_d(x1, -2.5, MPFR_RNDN);
            mpfr_set_d(y1, -1.5, MPFR_RNDN);
            mpfr_set_d(x2, 1.0, MPFR_RNDN);
            mpfr_set_d(y2, 1.5, MPFR_RNDN);
        }

        frame_count++;
        current_time = time(NULL);
        double elapsed = difftime(current_time, start_time);

        if (elapsed >= 1.0) {
            const char* zoom_mode = auto_zoom ? "Auto" : "Manual";
            double tx = mpfr_get_d(targetX, MPFR_RNDN);
            double ty = mpfr_get_d(targetY, MPFR_RNDN);
            printf("\rZoom: %.2e×   Mode: %s   Position: (%.15g, %.15g)   Iterations: %d   FPS: %.1f   ",
                   zoom_level_d, zoom_mode, tx, ty, max_iter, frame_count / elapsed);
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
    if (ref_real) {
        for (int i = 0; i < MAX_REFERENCE_ITERATIONS; i++) {
            mpfr_clear(ref_real[i]);
            mpfr_clear(ref_imag[i]);
            mpfr_clear(ref_dz_real[i]);
            mpfr_clear(ref_dz_imag[i]);
        }
        free(ref_real); free(ref_imag); free(ref_dz_real); free(ref_dz_imag);
    }
    munmap(fbp, screen_size);
    close(fb_fd);

    mpfr_clears(targetX, targetY, x1, y1, x2, y2, width, height, initial_width, temp, (mpfr_ptr)0);
    mpfr_clears(ref_center_re, ref_center_im, (mpfr_ptr)0);
    mpfr_clears(zoom_level, zoom_speed, min_width, manual_zoom_factor, (mpfr_ptr)0);

    tcsetattr(STDIN_FILENO, TCSANOW, &old_tio);
    return 0;
}
