#include "utils.h"

#define GRID_U_SAMPLES   80
#define GRID_V_SAMPLES   40
#define NUM_STREAMLINES  12
#define STREAMLINE_STEPS 200

#define R_MAJOR 2.0
#define r_MINOR 1.0

// Parametric torus
void torus_point(double u, double v, double *x, double *y, double *z)
{
  double cu = cos(u), su = sin(u);
  double cv = cos(v), sv = sin(v);
  double w  = R_MAJOR + r_MINOR * cv;
  *x = w * cu;
  *y = w * su;
  *z = r_MINOR * sv;
}

// One frame: sublevel + full streamlines
void draw_torus_sublevel_streamlines(double time)
{
  clear_framebuffer();

  // rotate torus in two planes
  double a = time * 0.3, b = time * 0.5;
  double ca = cos(a), sa = sin(a);
  double cb = cos(b), sb = sin(b);

  // compute animated threshold H(t) in [−r_MINOR, +r_MINOR]
  double H = r_MINOR * sin(time * 0.4);

  // draw sublevel {rotated_z <= H}
  for (int iu = 0; iu < GRID_U_SAMPLES; iu++) {
    double u = 2*PI * iu / GRID_U_SAMPLES;
    for (int iv = 0; iv < GRID_V_SAMPLES; iv++) {
      double v = 2*PI * iv / GRID_V_SAMPLES;
      double x,y,z;
      torus_point(u, v, &x, &y, &z);
      // rotate (x,y,z) → (x1,y1,z1)
      double y1 =  ca*y - sa*z;
      double z1 =  sa*y + ca*z;
      double x1 =  cb*x - sb*z1;
      double z2 =  sb*x + cb*z1;
      // if inside sublevel, paint a translucent pixel
      if (z2 <= H) {
        int px, py; double d;
        project_point(x1, y1, z2, &px, &py, &d);
        // pale yellow, alpha = 0.3
        draw_pixel(px, py, 255, 240, 128, 0.8);
      }
    }
  }

  // draw full streamlines of α = du/(2π) (horizontal loops, red)
  for (int k = 0; k < NUM_STREAMLINES; k++) {
    double v0 = 2*PI * k / NUM_STREAMLINES;
    int prev_px = -1, prev_py = -1;
    for (int i = 0; i <= STREAMLINE_STEPS; i++) {
      double u = 2*PI * i / STREAMLINE_STEPS + time * 0.1;
      double x,y,z;
      torus_point(u, v0, &x, &y, &z);
      double y1 =  ca*y - sa*z;
      double z1 =  sa*y + ca*z;
      double x1 =  cb*x - sb*z1;
      double z2 =  sb*x + cb*z1;
      int px, py; double d;
      project_point(x1, y1, z2, &px, &py, &d);
      rgb_color col = hsv_to_rgb(0.0, 0.8, 0.9); // red-ish
      if (prev_px>=0) {
        // draw a little segment by interpolating pixels
        draw_pixel(prev_px, prev_py, col.r, col.g, col.b, 1.0);
        draw_pixel(px,      py,      col.r, col.g, col.b, 1.0);
      }
      prev_px = px; prev_py = py;
    }
  }

  // 4) draw full streamlines of β = dv/(2π) (vertical loops, blue)
  for (int k = 0; k < NUM_STREAMLINES; k++) {
    double u0 = 2*PI * k / NUM_STREAMLINES;
    int prev_px = -1, prev_py = -1;
    for (int i = 0; i <= STREAMLINE_STEPS; i++) {
      double v = 2*PI * i / STREAMLINE_STEPS - time * 0.07;
      double x,y,z;
      torus_point(u0, v, &x, &y, &z);
      double y1 =  ca*y - sa*z;
      double z1 =  sa*y + ca*z;
      double x1 =  cb*x - sb*z1;
      double z2 =  sb*x + cb*z1;
      int px, py; double d;
      project_point(x1, y1, z2, &px, &py, &d);
      rgb_color col = hsv_to_rgb(220.0, 0.7, 0.9); // blue-ish
      if (prev_px>=0) {
        draw_pixel(prev_px, prev_py, col.r, col.g, col.b, 1.0);
        draw_pixel(px,      py,      col.r, col.g, col.b, 1.0);
      }
      prev_px = px; prev_py = py;
    }
  }
}

int main()
{
  if (!init_framebuffer()) return 1;
  double t = 0.0;
  while (1) {
    draw_torus_sublevel_streamlines(t);
    t += ROTATION_SPEED;
    usleep(16667);
  }
  cleanup_framebuffer();
  return 0;
}
