// dec_mesh_rotate.c
#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#define NX         40
#define NY         20
#define PI         3.14159265358979323846
#define R_MAJOR    2.0
#define r_MINOR    0.7
#define ROT_SPEED  0.02

#define MAX_EDGES  (NX*NY*3)
#define MAX_FACES  (NX*NY*2)

typedef struct { int v1, v2; } Edge;
typedef struct { int v[3], e[3], inc[3]; } Face;

// vertex‐to‐grid coords
static int    v_i[NX*NY], v_j[NX*NY];
static double VX[NX*NY], VY[NX*NY], VZ[NX*NY];

// mesh connectivity
static Edge   edges[MAX_EDGES];
static Face   faces[MAX_FACES];
static int    edge_count = 0, face_count = 0;

static void draw_line(int x0, int y0, int x1, int y1, uint8_t r, uint8_t g, uint8_t b) {
  int dx = x1 - x0, dy = y1 - y0;
  int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);
  if (steps == 0) {
    draw_pixel(x0, y0, r, g, b, 1.0);
    return;
  }
  double sx = dx / (double)steps, sy = dy / (double)steps;
  double x = x0, y = y0;
  for (int i = 0; i <= steps; i++) {
    draw_pixel((int)(x+0.5), (int)(y+0.5), r, g, b, 1.0);
    x += sx; y += sy;
  }
}

static void draw_circle(int xc, int yc, int rad, uint8_t r, uint8_t g, uint8_t b) {
  for (int dx = -rad; dx <= rad; dx++) {
    for (int dy = -rad; dy <= rad; dy++) {
      if (dx*dx + dy*dy <= rad*rad) {
        draw_pixel(xc+dx, yc+dy, r, g, b, 1.0);
      }
    }
  }
}

// rotate point (x,y,z) by yaw around Y then pitch around X
static inline void rotate3(double x, double y, double z,
                           double cosY, double sinY,
                           double cosX, double sinX,
                           double *xo, double *yo, double *zo) {
  // yaw:
  double x1 =  cosY*x + sinY*z;
  double z1 = -sinY*x + cosY*z;
  // pitch:
  double y1 =  cosX*y - sinX*z1;
  double z2 =  sinX*y + cosX*z1;
  *xo = x1; *yo = y1; *zo = z2;
}

// embed torus point in R^3
static void torus_point(double u, double v,
                        double *x, double *y, double *z) {
  double cu = cos(u), su = sin(u);
  double cv = cos(v), sv = sin(v);
  double w  = R_MAJOR + r_MINOR * cv;
  *x = w * cu;
  *y = w * su;
  *z = r_MINOR * sv;
}

int main() {
  if (!init_framebuffer()) return 1;

  // build vertices on torus
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      int id = i*NY + j;
      v_i[id] = i; v_j[id] = j;
      double u = 2*PI * i / NX;
      double v = 2*PI * j / NY;
      torus_point(u, v, &VX[id], &VY[id], &VZ[id]);
    }
  }

  // build faces and unique edges
  edge_count = face_count = 0;
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      int i2 = (i+1)%NX, j2 = (j+1)%NY;
      int v0 = i*NY + j,
          v1 = i2*NY + j,
          v2 = i*NY + j2,
          v3 = i2*NY + j2;
      int tris[2][3] = {{v0,v1,v2},{v1,v3,v2}};
      for (int t = 0; t < 2; t++) {
        Face *F = &faces[face_count++];
        for (int k = 0; k < 3; k++)
          F->v[k] = tris[t][k];
        for (int k = 0; k < 3; k++) {
          int A = F->v[k], B = F->v[(k+1)%3];
          int a = A<B?A:B, b = A<B?B:A;
          int eidx = -1;
          for (int e = 0; e < edge_count; e++)
            if (edges[e].v1==a && edges[e].v2==b) { eidx = e; break; }
          if (eidx<0) {
            eidx = edge_count;
            edges[edge_count].v1 = a;
            edges[edge_count].v2 = b;
            edge_count++;
          }
          F->e[k]   = eidx;
          F->inc[k] = (A==a && B==b)? +1 : -1;
        }
      }
    }
  }

  // allocate DEC arrays
  double *f0 = malloc(sizeof(double)*(NX*NY));
  double *a1 = malloc(sizeof(double)*edge_count);
  double *b2 = malloc(sizeof(double)*face_count);

  double t = 0.0;
  while (1) {
    // build 0-form f0
    for (int v = 0; v < NX*NY; v++) {
      double u = 2*PI * v_i[v]/NX;
      double w = 2*PI * v_j[v]/NY;
      f0[v] = sin(u + t) + cos(w + 1.3*t);
    }

    // d0 f0 -> a1 on edges (+ a little swirl on vertical edges)
    double maxA = 1e-9;
    for (int e = 0; e < edge_count; e++) {
      int i = edges[e].v1, j = edges[e].v2;
      double df = f0[j] - f0[i];
      // add small harmonic swirl on "vertical" edges:
      int di = (v_i[j] + NX - v_i[i])%NX;
      int dj = (v_j[j] + NY - v_j[i])%NY;
      double swirl = 0.0;
      if (di==0 && (dj==1||dj==NY-1)) {
        double u = 2*PI * v_i[i]/NX;
        double s = 0.3 * sin(u + 2.0*t);
        swirl = ((dj==1)? +s : -s);
      }
      a1[e] = df + swirl;
      maxA = fmax(maxA, fabs(a1[e]));
    }
    if (maxA<1e-6) maxA=1.0;

    // d1 a1 -> b2 on faces
    double maxB = 1e-9;
    for (int f = 0; f < face_count; f++) {
      double sum = 0;
      for (int k = 0; k < 3; k++)
        sum += faces[f].inc[k] * a1[faces[f].e[k]];
      b2[f] = sum;
      maxB = fmax(maxB, fabs(sum));
    }
    if (maxB<1e-6) maxB=1.0;

    // compute rotation angles
    double angleY = t * 0.3, angleX = t * 0.2;
    double cY = cos(angleY), sY = sin(angleY);
    double cX = cos(angleX), sX = sin(angleX);

    clear_framebuffer();

    // draw 2-form b2 on faces as colored dots
    for (int f = 0; f < face_count; f++) {
      int v0 = faces[f].v[0],
          v1 = faces[f].v[1],
          v2 = faces[f].v[2];
      double cx = (VX[v0]+VX[v1]+VX[v2])/3.0;
      double cy = (VY[v0]+VY[v1]+VY[v2])/3.0;
      double cz = (VZ[v0]+VZ[v1]+VZ[v2])/3.0;
      double xr, yr, zr;
      rotate3(cx,cy,cz, cY,sY, cX,sX, &xr,&yr,&zr);
      int px, py; double d;
      project_point(xr, yr, zr, &px, &py, &d);
      double c = b2[f]/maxB;
      if (c>1) c=1; 
      if (c<-1) c=-1;
      uint8_t R,G,B;
      if (c>=0) { R=0;    G=(uint8_t)(c*255); B=0; }
      else      { R=(uint8_t)(-c*255); G=0;    B=0; }
      draw_circle(px, py, 3, R, G, B);
    }

    // draw 1-form a1 on edges as colored arrows
    for (int e = 0; e < edge_count; e++) {
      int i = edges[e].v1, j = edges[e].v2;
      double x1 = VX[i], y1 = VY[i], z1 = VZ[i];
      double x2 = VX[j], y2 = VY[j], z2 = VZ[j];
      double xm = (x1+x2)*0.5, ym = (y1+y2)*0.5, zm = (z1+z2)*0.5;
      double dx = x2-x1, dy = y2-y1, dz = z2-z1;
      double L = sqrt(dx*dx+dy*dy+dz*dz);
      if (L<1e-6) continue;
      dx/=L; dy/=L; dz/=L;
      double xt = xm + dx*0.3*L,
             yt = ym + dy*0.3*L,
             zt = zm + dz*0.3*L;
      double xmr,ymr,zmr, xtr,ytr,ztr;
      rotate3(xm, ym, zm, cY,sY, cX,sX, &xmr,&ymr,&zmr);
      rotate3(xt, yt, zt, cY,sY, cX,sX, &xtr,&ytr,&ztr);
      int x0,y0,x1_,y1_; double d0,d1;
      project_point(xmr, ymr, zmr, &x0, &y0, &d0);
      project_point(xtr, ytr, ztr, &x1_,&y1_,&d1);
      double c = a1[e]/maxA;
      if (c>1) c=1;
      if (c<-1) c=-1;
      uint8_t R,G,B;
      if (c>=0) { R=(uint8_t)(c*255); G=0;            B=0; }
      else      { R=0;            G=0;            B=(uint8_t)(-c*255); }
      draw_line(x0,y0,x1_,y1_, R,G,B);
      // arrow head
      double ux=x1_-x0, uy=y1_-y0;
      double L2=sqrt(ux*ux+uy*uy);
      if (L2>0) {
        ux/=L2; uy/=L2;
        double pxp=-uy, pyp=ux;
        int hl=4;
        int hx1=(int)(x1_-ux*hl+pxp*hl+0.5),
            hy1=(int)(y1_-uy*hl+pyp*hl+0.5),
            hx2=(int)(x1_-ux*hl-pxp*hl+0.5),
            hy2=(int)(y1_-uy*hl-pyp*hl+0.5);
        draw_line(x1_,y1_, hx1,hy1, R,G,B);
        draw_line(x1_,y1_, hx2,hy2, R,G,B);
      }
    }

    // draw 0-form f0 at vertices as gray dots
    for (int v = 0; v < NX*NY; v++) {
      double xv,yv,zv;
      rotate3(VX[v], VY[v], VZ[v], cY,sY, cX,sX, &xv,&yv,&zv);
      int px, py; double d;
      project_point(xv, yv, zv, &px, &py, &d);
      double norm = (f0[v]+2.0)*0.25;
      if (norm<0) norm=0;
      if (norm>1) norm=1;
      uint8_t val = (uint8_t)(norm*255);
      draw_circle(px, py, 2, val, val, val);
    }

    t += ROT_SPEED;
    usleep(16667);
  }

  // never reached
  cleanup_framebuffer();
  free(f0); free(a1); free(b2);
  return 0;
}
