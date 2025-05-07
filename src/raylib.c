#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <time.h>
#include "raylib.h"

#define TARGET_RE -0.743643887037158704752191506114774
#define TARGET_IM 0.131825904205311970493132056385139

#define ZOOM_SPEED 0.97f
#define BASE_ITERATIONS 100
#define ITERATION_FACTOR 100
#define THREAD_COUNT 8
#define MAX_REFERENCE_ITERATIONS 10000
#define GLITCH_THRESHOLD 1e-12

#define AUTO_ZOOM_TIMEOUT 2
#define MANUAL_ZOOM_FACTOR 0.5
#define MOVEMENT_FACTOR 0.05

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720

double *ref_real;
double *ref_imag;
double *ref_dz_real;
double *ref_dz_imag;
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

static inline double calculate_smooth_value(double zr, double zi) {
    double mag_squared = zr*zr + zi*zi;
    return 1.0 - log(log(mag_squared) / 2.0) / LOG2;
}

static inline void calculate_reference_orbit(double center_re, double center_im, int max_iter) {
    printf("Calculating reference orbit for (%g, %g)...\n", center_re, center_im);

    if (!ref_real || !ref_imag || !ref_dz_real || !ref_dz_imag) {
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
    double dzr = 0.0, dzi = 0.0;

    int i;
    for (i = 0; i < max_iter; i++) {
        ref_real[i] = zr;
        ref_imag[i] = zi;
        ref_dz_real[i] = dzr;
        ref_dz_imag[i] = dzi;

        if (zr*zr + zi*zi > 4.0)
            break;

        double dz_temp = 2.0 * (zr * dzr - zi * dzi) + 1.0;
        dzi = 2.0 * (zr * dzi + zi * dzr);
        dzr = dz_temp;

        double temp = zr*zr - zi*zi + center_re;
        zi = 2.0 * zr * zi + center_im;
        zr = temp;
    }

    ref_iterations = i;
    printf("Reference orbit calculated with %d iterations\n", ref_iterations);
}

typedef struct {
    int width, height;
    double x1, y1, x2, y2;
    double center_re, center_im;
    int max_iter;
    int start_y, end_y;
    unsigned int *thread_buffer;
} ThreadData;

static inline void* render_slice(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    int width = data->width;
    double dx = (data->x2 - data->x1) / width;
    double dy = (data->y2 - data->y1) / data->height;

    for (int y = data->start_y; y < data->end_y; y++) {
        double cy = data->y1 + y * dy;
        int buffer_y = y - data->start_y;
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
                double zr = ref_zr + delta_zr;
                double zi = ref_zi + delta_zi;
                double mag_squared = zr*zr + zi*zi;
                if (mag_squared > 4.0) {
                    escaped = 1;
                    smooth_val = calculate_smooth_value(zr, zi);
                    break;
                }
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
                double temp = 2.0 * (ref_zr * delta_zr - ref_zi * delta_zi) +
                              delta_zr*delta_zr - delta_zi*delta_zi + delta_re;
                delta_zi = 2.0 * (ref_zr * delta_zi + ref_zi * delta_zr) +
                           2.0 * delta_zr * delta_zi + delta_im;
                delta_zr = temp;
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
    return NULL;
}

static inline void combine_thread_buffers(ThreadData* thread_data, int thread_count, int width, int height) {
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

static inline void update_texture_from_draw_buffer(Color *pixels, int width, int height) {
    for (int i = 0; i < width * height; i++) {
        unsigned int c = draw_buffer[i];
        pixels[i].r = (c >> 16) & 0xFF;
        pixels[i].g = (c >> 8) & 0xFF;
        pixels[i].b = c & 0xFF;
        pixels[i].a = 255;
    }
}

int main() {
    int width = WINDOW_WIDTH;
    int height = WINDOW_HEIGHT;

    InitWindow(width, height, "Mandelbrot Explorer (raylib)");
    SetTargetFPS(60);

    draw_buffer = malloc(width * height * sizeof(unsigned int));
    if (!draw_buffer) {
        perror("Error allocating drawing buffer");
        return 1;
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

    int slice_height = height / THREAD_COUNT;
    for (int i = 0; i < THREAD_COUNT; i++) {
        int end_y = (i == THREAD_COUNT - 1) ? height : (i + 1) * slice_height;
        int buffer_height = end_y - (i * slice_height);
        thread_data[i].thread_buffer = malloc(width * buffer_height * sizeof(unsigned int));
        if (!thread_data[i].thread_buffer) {
            perror("Error allocating thread buffer");
            goto end;
        }
    }

    Color *pixels = malloc(width * height * sizeof(Color));
    Texture2D texture = LoadTextureFromImage(GenImageColor(width, height, BLACK));

    int auto_zoom = 1;
    time_t last_input_time = time(NULL);

    printf("Mandelbrot Explorer with Navigation (raylib)\n");
    printf("Arrow keys: move | Z: zoom in | X: zoom out | C: toggle auto zoom | R: recenter | ESC: exit\n");
    printf("Auto-zoom resumes after 2 seconds of inactivity.\n");

    int running = 1;
    int frame_count = 0;
    time_t start_time = time(NULL);
    double initial_width = 3.5;
    int recalculate_reference = 0;

    while (!WindowShouldClose() && running) {
        double width_c = x2 - x1;
        double height_c = y2 - y1;

        // --- Input handling ---
        int moved = 0;
        double move_x = width_c * MOVEMENT_FACTOR;
        double move_y = height_c * MOVEMENT_FACTOR;

        if (IsKeyDown(KEY_RIGHT)) { targetX += move_x; moved = 1; }
        if (IsKeyDown(KEY_LEFT))  { targetX -= move_x; moved = 1; }
        if (IsKeyDown(KEY_UP))    { targetY -= move_y; moved = 1; }
        if (IsKeyDown(KEY_DOWN))  { targetY += move_y; moved = 1; }
        if (IsKeyPressed(KEY_Z))  { width_c *= MANUAL_ZOOM_FACTOR; height_c *= MANUAL_ZOOM_FACTOR; moved = 1; }
        if (IsKeyPressed(KEY_X))  { width_c /= MANUAL_ZOOM_FACTOR; height_c /= MANUAL_ZOOM_FACTOR; moved = 1; }
        if (IsKeyPressed(KEY_C))  { auto_zoom = !auto_zoom; }
        if (IsKeyPressed(KEY_R))  { targetX = TARGET_RE; targetY = TARGET_IM; width_c = initial_width; height_c = initial_width * height / width; moved = 1; }
        if (IsKeyPressed(KEY_ESCAPE)) { running = 0; break; }

        if (moved) {
            recalculate_reference = 1;
            last_input_time = time(NULL);
        }

        x1 = targetX - width_c / 2;
        x2 = targetX + width_c / 2;
        y1 = targetY - height_c / 2;
        y2 = targetY + height_c / 2;

        time_t current_time = time(NULL);
        if (!auto_zoom && difftime(current_time, last_input_time) >= AUTO_ZOOM_TIMEOUT) {
            auto_zoom = 1;
            printf("\rAuto-zoom resumed                                      \n");
        }

        if (recalculate_reference) {
            calculate_reference_orbit(targetX, targetY, MAX_REFERENCE_ITERATIONS);
            recalculate_reference = 0;
        }

        double zoom_level = initial_width / (x2 - x1);

        int max_iter = BASE_ITERATIONS + (int)(ITERATION_FACTOR * log10(zoom_level));
        if (max_iter < BASE_ITERATIONS) max_iter = BASE_ITERATIONS;
        if (max_iter > MAX_REFERENCE_ITERATIONS) max_iter = MAX_REFERENCE_ITERATIONS;

        for (int i = 0; i < THREAD_COUNT; i++) {
            thread_data[i].width = width;
            thread_data[i].height = height;
            thread_data[i].x1 = x1;
            thread_data[i].y1 = y1;
            thread_data[i].x2 = x2;
            thread_data[i].y2 = y2;
            thread_data[i].center_re = targetX;
            thread_data[i].center_im = targetY;
            thread_data[i].max_iter = max_iter;
            thread_data[i].start_y = i * slice_height;
            thread_data[i].end_y = (i == THREAD_COUNT - 1) ? height : (i + 1) * slice_height;
            pthread_create(&threads[i], NULL, render_slice, &thread_data[i]);
        }
        for (int i = 0; i < THREAD_COUNT; i++) {
            pthread_join(threads[i], NULL);
        }
        combine_thread_buffers(thread_data, THREAD_COUNT, width, height);

        update_texture_from_draw_buffer(pixels, width, height);
        UpdateTexture(texture, pixels);

        BeginDrawing();
        ClearBackground(BLACK);
        DrawTexture(texture, 0, 0, WHITE);

        // UI overlay
        double elapsed = difftime(current_time, start_time);
        frame_count++;
        if (elapsed >= 1.0) {
            start_time = current_time;
            frame_count = 0;
        }
        char info[256];
        snprintf(info, sizeof(info),
            "Zoom: %.2e×   Mode: %s   Position: (%.15g, %.15g)   Iterations: %d   FPS: %.1f",
            zoom_level, auto_zoom ? "Auto" : "Manual", targetX, targetY, max_iter, GetFPS());
        DrawText(info, 10, 10, 18, RAYWHITE);

        EndDrawing();

        if (auto_zoom && width_c > 1e-14) {
            double newWidth = width_c * ZOOM_SPEED;
            double newHeight = height_c * ZOOM_SPEED;
            x1 = targetX - newWidth / 2;
            x2 = targetX + newWidth / 2;
            y1 = targetY - newHeight / 2;
            y2 = targetY + newHeight / 2;
        } else if (width_c <= 1e-14) {
            x1 = -2.5;
            y1 = -1.5;
            x2 = 1.0;
            y2 = 1.5;
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
    free(pixels);
    UnloadTexture(texture);
    CloseWindow();
    return 0;
}
