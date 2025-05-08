#include <complex.h>
#include <stdio.h>
#include "utils.h"

// Complex matrix operations for U(n) representations
typedef struct {
    double complex data[2][2];
} complex_matrix2;

complex_matrix2 matrix_multiply(complex_matrix2 a, complex_matrix2 b) {
    complex_matrix2 result = {{{0}}};
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}

double complex u1_element(double t) {
    return cos(t) + I * sin(t);
}

// Create a U(2) matrix parameterized by angles
complex_matrix2 u2_element(double t, double s, double r) {
    complex_matrix2 matrix;
    matrix.data[0][0] = cos(t) * cexp(I * s);
    matrix.data[0][1] = sin(t) * cexp(I * r);
    matrix.data[1][0] = -sin(t) * cexp(-I * r);
    matrix.data[1][1] = cos(t) * cexp(-I * s);
    return matrix;
}

// Visualize Bott periodicity by showing loop spaces of U(1) and U(2)
void visualize_bott_periodicity(double time) {
    clear_framebuffer();
    
    double angle = time * 0.5;
    double x, y, z;
    int px, py;
    double depth;
    rgb_color color;
    
    int currentPeriod = ((int)(time / 10.0)) % 2;  // Toggle between 0 and 1 (U(1) and U(2))
    
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            // Map grid to parameters
            double u = 2.0 * PI * i / GRID_SIZE;
            double v = 2.0 * PI * j / GRID_SIZE;
            
            if (currentPeriod == 0) {  // Visualize π₁(U(1)) ≅ Z
                // In this case, we visualize loops in U(1) as circular paths
                // U(1) can be visualized as a circle in 2D
                double complex z1 = u1_element(u);
                
                x = creal(z1) * cos(angle);
                y = creal(z1) * sin(angle);
                z = cimag(z1);
                
                // Color based on winding number (integral part of u/(2π))
                double windingNumber = u / (2.0 * PI);
                color = hsv_to_rgb(fmod(windingNumber * 360.0, 360.0), 0.8, 0.9);
            } 
            else {  // Visualize π₃(U(2)) ≅ Z
                // We visualize loops in U(2)
                complex_matrix2 mat = u2_element(u, v, angle);
                
                double complex det = mat.data[0][0] * mat.data[1][1] - mat.data[0][1] * mat.data[1][0];
                
                x = creal(det) * cos(v);
                y = creal(det) * sin(v);
                z = cimag(det);
                
                // Color based on topological degree
                double phase = carg(det) / (2.0 * PI);
                color = hsv_to_rgb(fmod(phase * 360.0, 360.0), 0.8, 0.9);
            }
            
            // Apply transition effect between periods
            if (fmod(time, 10.0) < 5.0) {
                // Transitioning from period 0 to 1
                double transitionFactor = fmod(time, 5.0) / 5.0;
                x = x * (1.0 - transitionFactor) + (x * cos(v)) * transitionFactor;
                y = y * (1.0 - transitionFactor) + (y * sin(v)) * transitionFactor;
                z = z * (1.0 - transitionFactor) + (z + 0.2 * sin(u + v)) * transitionFactor;
            }
            
            project_point(x, y, z, &px, &py, &depth);
            
            double brightness = 0.5 + 0.5 * (depth + 1.5) / 3.0;
            draw_pixel(px, py, (uint8_t)(color.r * brightness), (uint8_t)(color.g * brightness), (uint8_t)(color.b * brightness), 1);
        }
    }
    
    printf("\r\033[2K");
    if (currentPeriod == 0) {
        printf("Visualizing pi_1(U(1)) = Z  [Period 0]");
    } else {
        printf("Visualizing pi_3(U(2)) = Z  [Period 1]");
    }
    fflush(stdout);
}

int main() {
    if (!init_framebuffer()) {
        fprintf(stderr, "Failed to initialize framebuffer\n");
        return 1;
    }
    
    printf("Bott Periodicity Visualization\n");
    printf("Press Ctrl+C to exit\n");
    
    double time = 0.0;
    while (1) {
        visualize_bott_periodicity(time);
        time += ROTATION_SPEED;
        usleep(16667);
    }
    
    cleanup_framebuffer();
    return 0;
}
