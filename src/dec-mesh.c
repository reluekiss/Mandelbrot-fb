// dec_heat_forms.c
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
#define ROT_SPEED  0.02  // rotation speed (rad/frame)

#define MAX_EDGES  (NX*NY*3)
#define MAX_FACES  (NX*NY*2)

typedef struct { int v1, v2; } Edge;
typedef struct { int v[3], e[3], inc[3]; } Face;

// vertex grid coords + embedded positions
static int    v_i[NX*NY], v_j[NX*NY];
static double VX[NX*NY], VY[NX*NY], VZ[NX*NY];

// mesh connectivity
static Edge   edges[MAX_EDGES];
static Face   faces[MAX_FACES];
static int    edge_count, face_count;

// simple line‐drawing
static void draw_line(int x0, int y0, int x1, int y1,
                      uint8_t r, uint8_t g, uint8_t b) {
  int dx = x1 - x0, dy = y1 - y0;
  int steps = abs(dx)>abs(dy)?abs(dx):abs(dy);
  if (!steps) {
    draw_pixel(x0, y0, r, g, b, 1.0);
    return;
  }
  double sx = dx/(double)steps, sy = dy/(double)steps;
  double x = x0, y = y0;
  for (int i = 0; i <= steps; i++) {
    draw_pixel((int)(x+0.5), (int)(y+0.5), r, g, b, 1.0);
    x += sx; y += sy;
  }
}

// filled circle for faces & vertices
static void draw_circle(int xc, int yc, int r,
                        uint8_t R, uint8_t G, uint8_t B) {
  for (int dx = -r; dx <= r; dx++)
    for (int dy = -r; dy <= r; dy++)
      if (dx*dx+dy*dy <= r*r)
        draw_pixel(xc+dx, yc+dy, R, G, B, 1.0);
}

// torus embedding (u,v) |-> (x,y,z)
static void torus_point(double u, double v,
                        double *x, double *y, double *z) {
  double cu = cos(u), su = sin(u),
         cv = cos(v), sv = sin(v);
  double w  = R_MAJOR + r_MINOR * cv;
  *x = w * cu; *y = w * su; *z = r_MINOR * sv;
}

// 3D rotate: yaw about Y then pitch about X
static inline void rotate3(double x,double y,double z,
                           double cy,double sy,
                           double cx,double sx,
                           double *xo,double *yo,double *zo){
  double x1 =  cy*x + sy*z;
  double z1 = -sy*x + cy*z;
  double y1 =  cx*y - sx*z1;
  double z2 =  sx*y + cx*z1;
  *xo = x1; *yo = y1; *zo = z2;
}

int main(){
  if(!init_framebuffer()) return 1;

  // build vertex grid & embed on torus
  for(int i=0;i<NX;i++)for(int j=0;j<NY;j++){
    int id = i*NY + j;
    v_i[id]=i; v_j[id]=j;
    double u=2*PI*i/NX, v=2*PI*j/NY;
    torus_point(u,v,&VX[id],&VY[id],&VZ[id]);
  }

  // build faces and unique edges
  edge_count=face_count=0;
  for(int i=0;i<NX;i++)for(int j=0;j<NY;j++){
    int i2=(i+1)%NX, j2=(j+1)%NY;
    int v0=i*NY+j, v1=i2*NY+j, v2=i*NY+j2, v3=i2*NY+j2;
    int tris[2][3]={{v0,v1,v2},{v1,v3,v2}};
    for(int t=0;t<2;t++){
      Face *F = &faces[face_count++];
      for(int k=0;k<3;k++) F->v[k]=tris[t][k];
      for(int k=0;k<3;k++){
        int A=F->v[k], B=F->v[(k+1)%3];
        int a=A<B?A:B, b=A<B?B:A, eidx=-1;
        for(int e=0;e<edge_count;e++){
          if(edges[e].v1==a && edges[e].v2==b){
            eidx=e; break;
          }
        }
        if(eidx<0){
          eidx=edge_count;
          edges[edge_count].v1=a;
          edges[edge_count].v2=b;
          edge_count++;
        }
        F->e[k]=eidx;
        F->inc[k]=(A==a && B==b)?+1:-1;
      }
    }
  }

  // allocate DEC arrays
  double *a1   = malloc(sizeof(double)*edge_count),
         *b2   = malloc(sizeof(double)*face_count),
         *v0   = malloc(sizeof(double)*(NX*NY)),
         *lap1 = malloc(sizeof(double)*edge_count);

  // init random 1-form on edges
  for(int e=0;e<edge_count;e++)
    a1[e] = 2.0*(rand()/(double)RAND_MAX) - 1.0;

  double t = 0.0;
  const double dt = ROT_SPEED; // simulation step

  while(1){
    // compute b2 = d1 a1
    double maxB=1e-9;
    for(int f=0;f<face_count;f++){
      double s=0;
      for(int k=0;k<3;k++)
        s += faces[f].inc[k] * a1[faces[f].e[k]];
      b2[f]=s;
      maxB = fmax(maxB,fabs(s));
    }
    if(maxB<1e-6) maxB=1.0;

    // compute δ1 a1 at vertices
    for(int v=0;v<NX*NY;v++) v0[v]=0;
    for(int e=0;e<edge_count;e++){
      int i=edges[e].v1, j=edges[e].v2;
      v0[i] -= a1[e];
      v0[j] += a1[e];
    }

    // build Laplacian on 1-forms
    double maxL=1e-9;
    for(int e=0;e<edge_count;e++){
      int i=edges[e].v1, j=edges[e].v2;
      double t1 = v0[j] - v0[i];
      double t2 = 0;
      // δ2 b2
      for(int f=0;f<face_count;f++){
        for(int k=0;k<3;k++){
          if(faces[f].e[k]==e)
            t2 += faces[f].inc[k] * b2[f];
        }
      }
      double L = t1 + t2;
      lap1[e]=L;
      maxL = fmax(maxL,fabs(L));
    }
    if(maxL<1e-6) maxL=1.0;

    // Euler step
    for(int e=0;e<edge_count;e++)
      a1[e] -= dt * lap1[e] / maxL;

    // compute rotation
    double angY = t*0.3, angX = t*0.2;
    double cY=cos(angY), sY=sin(angY),
           cX=cos(angX), sX=sin(angX);

    clear_framebuffer();

    // draw b2 on faces
    for(int f=0;f<face_count;f++){
      int v0i=faces[f].v[0], v1i=faces[f].v[1],
          v2i=faces[f].v[2];
      double cx=(VX[v0i]+VX[v1i]+VX[v2i])/3.0,
             cy=(VY[v0i]+VY[v1i]+VY[v2i])/3.0,
             cz=(VZ[v0i]+VZ[v1i]+VZ[v2i])/3.0;
      double xr,yr,zr;
      rotate3(cx,cy,cz, cY,sY, cX,sX, &xr,&yr,&zr);
      int px,py; double d;
      project_point(xr,yr,zr,&px,&py,&d);
      double c = b2[f]/maxB;
      if(c>1) c=1; if(c<-1) c=-1;
      uint8_t R,G,B;
      if(c>=0){ R=0; G=(uint8_t)(c*255); B=0;}
      else    { R=(uint8_t)(-c*255); G=0; B=0;}
      draw_circle(px,py,3,R,G,B);
    }

    // draw a1 on edges
    for(int e=0;e<edge_count;e++){
      int i=edges[e].v1, j=edges[e].v2;
      double x1=VX[i],y1=VY[i],z1=VZ[i],
             x2=VX[j],y2=VY[j],z2=VZ[j];
      double xm=(x1+x2)*0.5, ym=(y1+y2)*0.5,
             zm=(z1+z2)*0.5;
      double dx=x2-x1, dy=y2-y1, dz=z2-z1;
      double L3=sqrt(dx*dx+dy*dy+dz*dz);
      if(L3<1e-6) continue;
      dx/=L3; dy/=L3; dz/=L3;
      double xt=xm+dx*0.3*L3,
             yt=ym+dy*0.3*L3,
             zt=zm+dz*0.3*L3;
      double xr0,yr0,zr0, xr1,yr1,zr1;
      rotate3(xm,ym,zm, cY,sY, cX,sX,&xr0,&yr0,&zr0);
      rotate3(xt,yt,zt, cY,sY, cX,sX,&xr1,&yr1,&zr1);
      int x0,y0,x1_,y1_; double d0,d1;
      project_point(xr0,yr0,zr0,&x0,&y0,&d0);
      project_point(xr1,yr1,zr1,&x1_,&y1_,&d1);
      double c = a1[e];
      if(c>1) c=1; if(c<-1)c=-1;
      uint8_t R,G,B;
      if(c>=0){ R=(uint8_t)(c*255); G=0; B=0;}
      else    { R=0; G=0; B=(uint8_t)(-c*255);}
      draw_line(x0,y0,x1_,y1_, R,G,B);
      // arrowhead
      double ux=x1_-x0, uy=y1_-y0;
      double l2=sqrt(ux*ux+uy*uy);
      if(l2>0){
        ux/=l2; uy/=l2;
        double pxp=-uy, pyp=ux; int hl=4;
        int hx1=(int)(x1_-ux*hl+pxp*hl+0.5),
            hy1=(int)(y1_-uy*hl+pyp*hl+0.5);
        int hx2=(int)(x1_-ux*hl-pxp*hl+0.5),
            hy2=(int)(y1_-uy*hl-pyp*hl+0.5);
        draw_line(x1_,y1_,hx1,hy1,R,G,B);
        draw_line(x1_,y1_,hx2,hy2,R,G,B);
      }
    }

    // advance time & frame
    t += ROT_SPEED;
    usleep(16667);
  }

  // never reached
  cleanup_framebuffer();
  free(a1); free(b2); free(v0); free(lap1);
  return 0;
}
