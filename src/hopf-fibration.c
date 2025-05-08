#include "utils.h"

#define FIBER_SAMPLES 100
#define SPHERE_SAMPLES 200

// Hopf map from S^3 to S^2
void hopf_map(double x1, double x2, double x3, double x4, double *y1, double *y2, double *y3) {
  *y1 = 2.0 * (x1 * x3 + x2 * x4);
  *y2 = 2.0 * (x2 * x3 - x1 * x4);
  *y3 = x1 * x1 + x2 * x2 - x3 * x3 - x4 * x4;
}

// Generate a point on S^3 that maps to a specific point on S^2 with a phase angle
void inverse_hopf_map(double y1, double y2, double y3, double t, double *x1, double *x2, double *x3, double *x4) {
  double s2 = 1.0 - y3;
  double denom = sqrt(2.0 * s2);
  
  if (fabs(denom) < 1e-10) {
    // Handle the special case where y3 is close to 1
    *x1 = cos(t);
    *x2 = sin(t);
    *x3 = 0;
    *x4 = 0;
  } else {
    *x1 = sqrt((1.0 + y3) / 2.0);
    *x2 = y1 / denom;
    *x3 = y2 / denom;
    *x4 = 0.0;
    
    // Rotate by angle t along the fiber
    double c = cos(t);
    double s = sin(t);
    double nx1 = c * (*x1) - s * (*x4);
    double nx4 = s * (*x1) + c * (*x4);
    double nx2 = c * (*x2) - s * (*x3);
    double nx3 = s * (*x2) + c * (*x3);
    
    *x1 = nx1;
    *x2 = nx2;
    *x3 = nx3;
    *x4 = nx4;
  }
}

void rotate_4d(double *x1, double *x2, double *x3, double *x4, double angle1, double angle2, double angle3) {
  // Rotate in the x1-x2 plane
  double c1 = cos(angle1);
  double s1 = sin(angle1);
  double nx1 = c1 * (*x1) - s1 * (*x2);
  double nx2 = s1 * (*x1) + c1 * (*x2);
  
  // Rotate in the x3-x4 plane
  double c2 = cos(angle2);
  double s2 = sin(angle2);
  double nx3 = c2 * (*x3) - s2 * (*x4);
  double nx4 = s2 * (*x3) + c2 * (*x4);
  
  // Rotate in the x1-x3 plane
  double c3 = cos(angle3);
  double s3 = sin(angle3);
  *x1 = c3 * nx1 - s3 * nx3;
  *x3 = s3 * nx1 + c3 * nx3;
  *x2 = nx2;
  *x4 = nx4;
}

void draw_hopf_fibration(double time) {
  clear_framebuffer();
  
  for (int i = 0; i < SPHERE_SAMPLES; i++) {
    double theta = PI * (1.0 + sqrt(5.0)) * i; // Golden angle increment
    double z = 2.0 * ((double)i / SPHERE_SAMPLES) - 1.0;
    double radius = sqrt(1.0 - z * z);
    double y1 = radius * cos(theta);
    double y2 = radius * sin(theta);
    double y3 = z;
    
    for (int j = 0; j < FIBER_SAMPLES; j++) {
      double t = 2.0 * PI * j / FIBER_SAMPLES;
      double x1, x2, x3, x4;
      
      inverse_hopf_map(y1, y2, y3, t + time, &x1, &x2, &x3, &x4);
      
      rotate_4d(&x1, &x2, &x3, &x4, time * 0.5, time * 0.3, time * 0.2);
      
      int px, py;
      double depth = x4;
      project_point(x1, x2, x3, &px, &py, &depth);
      
      double hue = 360.0 * (t / (2.0 * PI));
      double saturation = 0.8 + 0.2 * sin(time + 3.0 * t);
      double value = 0.7 + 0.3 * (1.0 + y3) / 2.0;
      
      rgb_color color = hsv_to_rgb(hue, saturation, value);
      
      draw_pixel(px, py, (uint8_t)color.r, (uint8_t)color.g, (uint8_t)color.b, 1);
    }
  }
}

int main() {
  if (!init_framebuffer()) {
    fprintf(stderr, "Failed to initialize framebuffer\n");
    return 1;
  }
  
  printf("Hopf Fibration Visualization\n");
  printf("Press Ctrl+C to exit\n");
  
  double time = 0.0;
  while (1) {
    draw_hopf_fibration(time);
    time += ROTATION_SPEED;
    usleep(16667);
  }
  
  cleanup_framebuffer();
  return 0;
}
