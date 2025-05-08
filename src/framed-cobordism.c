#include <stdlib.h>
#include "utils.h"

#define SAMPLES_U 80
#define SAMPLES_V 30
#define VECTOR_SCALE 0.15

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

// Function to compute surface points and normal vectors for cobordism
void compute_cobordism_point(double u, double v, double time, 
                            vector3 *point, vector3 *normal1, vector3 *normal2) {
    // Parameter u goes around the manifold (0 to 2π)
    // Parameter v goes between the two boundary manifolds (0 to 1)
    
    double theta = u * 2.0 * PI;
    
    // Animate the cobordism with time
    double anim_phase = fmod(time * 0.2, 2.0 * PI);
    
    // First boundary manifold (v = 0): a circle in the xy-plane
    vector3 start_point = {cos(theta), sin(theta), -1.0};
    
    // Second boundary manifold (v = 1): a more complex curve (trefoil knot)
    double a = 0.5;
    double b = 0.1;
    double c = 0.25;
    vector3 end_point = {
        (a + b * cos(3 * theta)) * cos(2 * theta),
        (a + b * cos(3 * theta)) * sin(2 * theta),
        1.0 + c * sin(3 * theta)
    };
    
    // Interpolate between manifolds to create cobordism
    double tween = v;
    
    // Apply some animation to the interpolation
    tween = tween + 0.1 * sin(anim_phase + theta * 2.0) * (1.0 - tween) * tween;
    
    // Compute point on the cobordism
    point->x = start_point.x * (1.0 - tween) + end_point.x * tween;
    point->y = start_point.y * (1.0 - tween) + end_point.y * tween;
    point->z = start_point.z * (1.0 - tween) + end_point.z * tween;
    
    // Compute tangent vectors (derivatives)
    double delta = 0.01;
    
    // Tangent in u direction
    double theta_plus = (u + delta) * 2.0 * PI;
    vector3 start_point_u_plus = {cos(theta_plus), sin(theta_plus), -1.0};
    vector3 end_point_u_plus = {
        (a + b * cos(3 * theta_plus)) * cos(2 * theta_plus),
        (a + b * cos(3 * theta_plus)) * sin(2 * theta_plus),
        1.0 + c * sin(3 * theta_plus)
    };
    
    vector3 point_u_plus;
    point_u_plus.x = start_point_u_plus.x * (1.0 - tween) + end_point_u_plus.x * tween;
    point_u_plus.y = start_point_u_plus.y * (1.0 - tween) + end_point_u_plus.y * tween;
    point_u_plus.z = start_point_u_plus.z * (1.0 - tween) + end_point_u_plus.z * tween;
    
    // Tangent in v direction
    double tween_plus = v + delta;
    tween_plus = tween_plus + 0.1 * sin(anim_phase + theta * 2.0) * (1.0 - tween_plus) * tween_plus;
    
    vector3 point_v_plus;
    point_v_plus.x = start_point.x * (1.0 - tween_plus) + end_point.x * tween_plus;
    point_v_plus.y = start_point.y * (1.0 - tween_plus) + end_point.y * tween_plus;
    point_v_plus.z = start_point.z * (1.0 - tween_plus) + end_point.z * tween_plus;
    
    // Compute tangent vectors
    vector3 tangent_u;
    tangent_u.x = (point_u_plus.x - point->x) / delta;
    tangent_u.y = (point_u_plus.y - point->y) / delta;
    tangent_u.z = (point_u_plus.z - point->z) / delta;
    
    vector3 tangent_v;
    tangent_v.x = (point_v_plus.x - point->x) / delta;
    tangent_v.y = (point_v_plus.y - point->y) / delta;
    tangent_v.z = (point_v_plus.z - point->z) / delta;
    
    // Normal vector is cross product of tangents
    *normal1 = vector_normalize(vector_cross(tangent_u, tangent_v));
    
    // Second normal vector - create a framing by rotating first normal
    // This simulates a nontrivial framing that varies along the cobordism
    double rotation_angle = time * 0.5 + theta + v * PI;
    double nx = normal1->x;
    double ny = normal1->y;
    double nz = normal1->z;
    
    // Create a perpendicular vector to normal1
    vector3 perpendicular;
    if (fabs(nx) < fabs(ny) && fabs(nx) < fabs(nz)) {
        perpendicular.x = 0;
        perpendicular.y = -nz;
        perpendicular.z = ny;
    } else if (fabs(ny) < fabs(nz)) {
        perpendicular.x = -nz;
        perpendicular.y = 0;
        perpendicular.z = nx;
    } else {
        perpendicular.x = -ny;
        perpendicular.y = nx;
        perpendicular.z = 0;
    }
    
    perpendicular = vector_normalize(perpendicular);
    
    // Rotate perpendicular vector around normal1 to create second normal
    c = cos(rotation_angle);
    double s = sin(rotation_angle);
    vector3 temp = vector_cross(*normal1, perpendicular);
    
    normal2->x = perpendicular.x * c + temp.x * s;
    normal2->y = perpendicular.y * c + temp.y * s;
    normal2->z = perpendicular.z * c + temp.z * s;
    
    *normal2 = vector_normalize(*normal2);
}

void rotate_view(vector3 *p, double view_angle_x, double view_angle_y, double view_angle_z) {
    double temp_y = p->y * cos(view_angle_x) - p->z * sin(view_angle_x);
    double temp_z = p->y * sin(view_angle_x) + p->z * cos(view_angle_x);
    p->y = temp_y;
    p->z = temp_z;
    
    double temp_x = p->x * cos(view_angle_y) + p->z * sin(view_angle_y);
    temp_z = -p->x * sin(view_angle_y) + p->z * cos(view_angle_y);
    p->x = temp_x;
    p->z = temp_z;
    
    temp_x = p->x * cos(view_angle_z) - p->y * sin(view_angle_z);
    temp_y = p->x * sin(view_angle_z) + p->y * cos(view_angle_z);
    p->x = temp_x;
    p->y = temp_y;
}

void visualize_framed_cobordism(double time) {
    clear_framebuffer();
    
    double view_angle_x = time * 0.1;
    double view_angle_y = sin(time * 0.05) * 0.5;
    double view_angle_z = 0.3;
    
    
    for (int i = 0; i < SAMPLES_U; i++) {
        for (int j = 0; j <= SAMPLES_V; j++) {
            double u = (double)i / SAMPLES_U;
            double v = (double)j / SAMPLES_V;
            
            vector3 point, normal1, normal2;
            compute_cobordism_point(u, v, time, &point, &normal1, &normal2);
            
            rotate_view(&point, view_angle_x, view_angle_y, view_angle_z);
            rotate_view(&normal1, view_angle_x, view_angle_y, view_angle_z);
            rotate_view(&normal2, view_angle_x, view_angle_y, view_angle_z);
            
            int px, py;
            double depth;
            project_point(point.x, point.y, point.z, &px, &py, &depth);
            
            rgb_color color = hsv_to_rgb(fmod(u * 360.0, 360.0), 0.7 + 0.3 * v, 0.9);
            double brightness = 0.3 + 0.7 * (depth + 2.0) / 4.0;
            
            draw_pixel(px, py, color.r * brightness, color.g * brightness, color.b * brightness, 1.0);
            
            vector3 normal1_end = vector_add(point, vector_scale(normal1, VECTOR_SCALE));
            vector3 normal2_end = vector_add(point, vector_scale(normal2, VECTOR_SCALE));
            
            int nx1, ny1, nx2, ny2;
            double ndepth1, ndepth2;
            project_point(normal1_end.x, normal1_end.y, normal1_end.z, &nx1, &ny1, &ndepth1);
            project_point(normal2_end.x, normal2_end.y, normal2_end.z, &nx2, &ny2, &ndepth2);
            
            if ((i % 5 == 0) && (j % 2 == 0)) {
                // First normal (red)
                draw_line(px, py, nx1, ny1, 255, 50, 50);
                
                // Second normal (blue) - showing the framing
                draw_line(px, py, nx2, ny2, 50, 50, 255);
            }
            
            // Connect to adjacent points to show the surface
            if (i > 0 && j > 0) {
                // Get the previous points
                vector3 prev_point_u, prev_normal1_u, prev_normal2_u;
                compute_cobordism_point(u - 1.0/SAMPLES_U, v, time, &prev_point_u, &prev_normal1_u, &prev_normal2_u);
                
                vector3 prev_point_v, prev_normal1_v, prev_normal2_v;
                compute_cobordism_point(u, v - 1.0/SAMPLES_V, time, &prev_point_v, &prev_normal1_v, &prev_normal2_v);
                
                rotate_view(&prev_point_u, view_angle_x, view_angle_y, view_angle_z);
                rotate_view(&prev_point_v, view_angle_x, view_angle_y, view_angle_z);
                
                int prev_px_u, prev_py_u, prev_px_v, prev_py_v;
                double prev_depth_u, prev_depth_v;
                project_point(prev_point_u.x, prev_point_u.y, prev_point_u.z, &prev_px_u, &prev_py_u, &prev_depth_u);
                project_point(prev_point_v.x, prev_point_v.y, prev_point_v.z, &prev_px_v, &prev_py_v, &prev_depth_v);
                
                // Draw lines to create surface grid
                if (i % 2 == 0)
                    draw_line(px, py, prev_px_u, prev_py_u, 100, 100, 100);
                if (j % 2 == 0)
                    draw_line(px, py, prev_px_v, prev_py_v, 100, 100, 100);
            }
        }
    }
}

int main() {
    if (!init_framebuffer()) {
        fprintf(stderr, "Failed to initialize framebuffer\n");
        return 1;
    }
    
    printf("Framed Cobordism: Manifolds with compatible framings\n");
    printf("Press Ctrl+C to exit\n");
    
    double time = 0.0;
    while (1) {
        visualize_framed_cobordism(time);
        time += ROTATION_SPEED;
        usleep(16667);
    }
    
    cleanup_framebuffer();
    return 0;
}
