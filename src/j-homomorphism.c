#include <stdlib.h>
#include "utils.h"

#define SAMPLES 200

typedef struct {
    double m[3][3];
} matrix3;

matrix3 so3_from_axis_angle(vector3 axis, double angle) {
    matrix3 r;
    double c = cos(angle);
    double s = sin(angle);
    double t = 1 - c;
    
    axis = vector_normalize(axis);
    
    r.m[0][0] = c + axis.x * axis.x * t;
    r.m[1][1] = c + axis.y * axis.y * t;
    r.m[2][2] = c + axis.z * axis.z * t;
    
    double tmp1 = axis.x * axis.y * t;
    double tmp2 = axis.z * s;
    r.m[0][1] = tmp1 - tmp2;
    r.m[1][0] = tmp1 + tmp2;
    
    tmp1 = axis.x * axis.z * t;
    tmp2 = axis.y * s;
    r.m[0][2] = tmp1 + tmp2;
    r.m[2][0] = tmp1 - tmp2;
    
    tmp1 = axis.y * axis.z * t;
    tmp2 = axis.x * s;
    r.m[1][2] = tmp1 - tmp2;
    r.m[2][1] = tmp1 + tmp2;
    
    return r;
}

vector3 matrix_apply(matrix3 m, vector3 v) {
    vector3 result;
    result.x = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z;
    result.y = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z;
    result.z = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z;
    return result;
}

vector3 hopf_map(double x1, double x2, double x3, double x4) {
    vector3 result;
    result.x = 2.0 * (x1 * x3 + x2 * x4);
    result.y = 2.0 * (x2 * x3 - x1 * x4);
    result.z = x1 * x1 + x2 * x2 - x3 * x3 - x4 * x4;
    return result;
}

void point_on_s3(double u, double v, double *x1, double *x2, double *x3, double *x4) {
    double theta1 = u * 2.0 * PI;
    double theta2 = v * 2.0 * PI;
    double phi = u * PI;
    
    *x1 = cos(phi) * cos(theta1);
    *x2 = cos(phi) * sin(theta1);
    *x3 = sin(phi) * cos(theta2);
    *x4 = sin(phi) * sin(theta2);
}

void draw_line(int x1, int y1, int x2, int y2, uint8_t r, uint8_t g, uint8_t b) {
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;
    
    while (1) {
        draw_pixel(x1, y1, r, g, b, 1.0);
        
        if (x1 == x2 && y1 == y2) break;
        
        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x1 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y1 += sy;
        }
    }
}

void visualize_j_homomorphism(double time) {
    clear_framebuffer();
    
    // Visualize the mapping from SO(3) to homotopy groups of spheres
    // We'll show how rotations in SO(3) map to elements in π_{n+3}(S^n)
    
    // First, draw a grid on S^3 (which we'll project to 3D)
    for (int i = 0; i < SAMPLES; i++) {
        for (int j = 0; j < SAMPLES; j += 5) {  // Skip some samples for clarity
            double u = (double)i / SAMPLES;
            double v = (double)j / SAMPLES;
            
            // Get point on S^3
            double x1, x2, x3, x4;
            point_on_s3(u, v, &x1, &x2, &x3, &x4);
            
            // Create a rotation (element of SO(3)) based on the point
            vector3 axis = {x1, x2, x3};
            double angle = x4 * 2.0 * PI + time;
            matrix3 rotation = so3_from_axis_angle(axis, angle);
            
            // Sample points around this rotation to show its action
            for (int k = 0; k < 12; k++) {
                double theta = k * 2.0 * PI / 12.0;
                
                vector3 base_vector = {cos(theta), sin(theta), 0.0};
                vector3 rotated = matrix_apply(rotation, base_vector);
                
                // The J-homomorphism sends this rotation to a map S^3 -> S^2
                // We visualize this by showing how the rotation transforms a sphere
                
                int px, py;
                double depth;
                
                double disp_x = rotated.x * 0.5 + x1 * 0.5;
                double disp_y = rotated.y * 0.5 + x2 * 0.5;
                double disp_z = rotated.z * 0.5 + x3 * 0.5;
                
                project_point(disp_x, disp_y, disp_z, &px, &py, &depth);
                
                // Calculate hue based on the "winding" of the J-homomorphism
                double winding = fmod(u * v * 8.0 + time * 0.1, 1.0);
                rgb_color color = hsv_to_rgb(winding * 360.0, 0.8, 0.9);
                
                draw_pixel(px, py, color.r, color.g, color.b, 0.5);
                
                // Draw connection to next point
                if (k > 0) {
                    int prev_px, prev_py;
                    double prev_depth;
                    double prev_theta = (k-1) * 2.0 * PI / 12.0;
                    vector3 prev_base = {cos(prev_theta), sin(prev_theta), 0.0};
                    vector3 prev_rotated = matrix_apply(rotation, prev_base);
                    
                    double prev_disp_x = prev_rotated.x * 0.5 + x1 * 0.5;
                    double prev_disp_y = prev_rotated.y * 0.5 + x2 * 0.5;
                    double prev_disp_z = prev_rotated.z * 0.5 + x3 * 0.5;
                    
                    project_point(prev_disp_x, prev_disp_y, prev_disp_z, &prev_px, &prev_py, &prev_depth);
                    
                    draw_line(prev_px, prev_py, px, py, color.r, color.g, color.b);
                }
            }
        }
    }
}

int main() {
    if (!init_framebuffer()) {
        fprintf(stderr, "Failed to initialize framebuffer\n");
        return 1;
    }
    
    printf("J-homomorphism: Mapping from SO(3) to homotopy groups of spheres\n");
    printf("Press Ctrl+C to exit\n");
    
    double time = 0.0;
    while (1) {
        visualize_j_homomorphism(time);
        time += ROTATION_SPEED;
        usleep(16667);
    }
    
    cleanup_framebuffer();
    return 0;
}
