#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <linux/fb.h>
#include <math.h>
#include <string.h>

struct fb_var_screeninfo vinfo;
struct fb_fix_screeninfo finfo;
uint8_t *fbp = NULL;
int fbfd = -1;
long screensize = 0;

#define PI 3.14159265358979323846
#define WIDTH 1024
#define HEIGHT 768
#define GRID_SIZE 200
#define ROTATION_SPEED 0.01

typedef struct {
  double r;
  double g;
  double b;
} rgb_color;

rgb_color hsv_to_rgb(double h, double s, double v) {
  rgb_color rgb;
  double c = v * s;
  double x = c * (1 - fabs(fmod(h / 60.0, 2) - 1));
  double m = v - c;

  if (h < 60) {
    rgb.r = c; rgb.g = x; rgb.b = 0;
  } else if (h < 120) {
    rgb.r = x; rgb.g = c; rgb.b = 0;
  } else if (h < 180) {
    rgb.r = 0; rgb.g = c; rgb.b = x;
  } else if (h < 240) {
    rgb.r = 0; rgb.g = x; rgb.b = c;
  } else if (h < 300) {
    rgb.r = x; rgb.g = 0; rgb.b = c;
  } else {
    rgb.r = c; rgb.g = 0; rgb.b = x;
  }

  rgb.r = (rgb.r + m) * 255;
  rgb.g = (rgb.g + m) * 255;
  rgb.b = (rgb.b + m) * 255;
  
  return rgb;
}

int init_framebuffer() {
  fbfd = open("/dev/fb0", O_RDWR);
  if (fbfd == -1) {
    perror("Error opening framebuffer device");
    return 0;
  }

  if (ioctl(fbfd, FBIOGET_FSCREENINFO, &finfo) == -1) {
    perror("Error reading fixed framebuffer information");
    close(fbfd);
    return 0;
  }

  if (ioctl(fbfd, FBIOGET_VSCREENINFO, &vinfo) == -1) {
    perror("Error reading variable framebuffer information");
    close(fbfd);
    return 0;
  }

  screensize = vinfo.xres * vinfo.yres * vinfo.bits_per_pixel / 8;
  fbp = (uint8_t *)mmap(0, screensize, PROT_READ | PROT_WRITE, MAP_SHARED, fbfd, 0);
  
  if ((long)fbp == -1) {
    perror("Error mapping framebuffer to memory");
    close(fbfd);
    return 0;
  }

  return 1;
}

void cleanup_framebuffer() {
  if (fbp != NULL) {
    munmap(fbp, screensize);
  }
  if (fbfd != -1) {
    close(fbfd);
  }
}

void draw_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, double alpha) {
    if (x < 0 || x >= vinfo.xres || y < 0 || y >= vinfo.yres) {
        return;
    }

    long location = (x + vinfo.xoffset) * (vinfo.bits_per_pixel / 8) +
                    (y + vinfo.yoffset) * finfo.line_length;

    if (vinfo.bits_per_pixel == 32) {
        uint8_t existing_b = *(fbp + location);
        uint8_t existing_g = *(fbp + location + 1);
        uint8_t existing_r = *(fbp + location + 2);
        
        *(fbp + location) = (uint8_t)(b * alpha + existing_b * (1 - alpha));
        *(fbp + location + 1) = (uint8_t)(g * alpha + existing_g * (1 - alpha));
        *(fbp + location + 2) = (uint8_t)(r * alpha + existing_r * (1 - alpha));
        *(fbp + location + 3) = 0;
    } else if (vinfo.bits_per_pixel == 16) {
        uint16_t pixel = ((r >> 3) << 11) | ((g >> 2) << 5) | (b >> 3);
        *((uint16_t *)(fbp + location)) = pixel;
    }
}

void clear_framebuffer() {
  memset(fbp, 0, screensize);
}

void project_point(double x, double y, double z, int *px, int *py, double *depth) {
  double scale = 300.0;
  double perspective = 1.0 / (2.0 - z);
  
  *px = vinfo.xres / 2 + (int)(x * scale * perspective);
  *py = vinfo.yres / 2 - (int)(y * scale * perspective);
  *depth = z;
}

typedef struct {
    double x;
    double y;
    double z;
} vector3;

// Vector operations
vector3 vector_cross(vector3 a, vector3 b) {
    vector3 result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

double vector_dot(vector3 a, vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vector3 vector_normalize(vector3 v) {
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len < 0.0001) {
        vector3 zero = {0, 0, 0};
        return zero;
    }
    
    vector3 result;
    result.x = v.x / len;
    result.y = v.y / len;
    result.z = v.z / len;
    return result;
}

vector3 vector_scale(vector3 v, double scale) {
    vector3 result;
    result.x = v.x * scale;
    result.y = v.y * scale;
    result.z = v.z * scale;
    return result;
}

vector3 vector_add(vector3 a, vector3 b) {
    vector3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

#endif // UTILS_H
