
#include <stdio.h>
#include <math.h>




extern int debug;

// The following class implements the EM convergence acceleration algorithm using Newton's
// iteration.

class EM_Newton
{
  double fx, fy; // next value for x,y;

  double *a, *b, *c; // a,b,c matrices as used in the EM
  int commark; // num of common typed markers

public:

  void init(double *a1, double *b1, double *c1, int commark1)
  {
    a = a1;
    b = b1;
    c = c1;
    commark = commark1;
  }

  // Performs one regular step of the EM algorithm
  
  void EM_getnext(double &x2, double &y2, double x, double y)
  {
    int i;
    double p0_next = 0, p1_next = 0, p2_next = 0;

    double p0 = x, p1 = y, p2 = 1-p0-p1;

    for (i = 1; i <= commark; i++) {
      p0_next += a[i] * p0 / (a[i] * p0 + b[i] * p1 + c[i] * p2);
      p1_next += b[i] * p1 / (a[i] * p0 + b[i] * p1 + c[i] * p2);
    }

    p0_next = p0_next/commark;
    p1_next = p1_next/commark;
    p2_next = 1.0 - p0_next-p1_next;

    x2 = p0_next;
    y2 = p1_next;
  }

  double diff(double x2, double y2, double x, double y)
  {
    return fabs(x2-x)+fabs(y2-y);
  }


  // The mapping function ensures the constraints for p0 and p1
  // it's not used right now

  void map(double &p0, double &p1)
  {
    //printf("map %f %f   ->  ",p0,p1);

    //p0 = (1+sin(p0))/2;
    //p1 = (1+sin(p1))/2;
    //p1 = p1*(1-p0);

    //p0 = (1+atan(p0)/1.570797)/2;
    //p1 = (1+atan(p1)/1.570797)/2;

    //p0 = (p0+2-sqrt(p0*p0+4))/(2*p0);
    //p1 = (p1+2-sqrt(p1*p1+4))/(2*p1);

    /*if(p0<2) p0 = 2+2-p0;
    if(p1<2) p1 = 2+2-p1;

    p0 = (p0-sqrt(p0*p0-4))/2;
    p1 = (p1-sqrt(p1*p1-4))/2;*/

    // printf("%f %f\n",p0,p1);
    //printf("p0,p1=%f %f\n",p0,p1);
  }

  // This is the function that we seek the zero for by Newton-Rhapson

  void optF(double &nX, double &nY, double p0, double p1)
  {
    // p0n <- 1-sum((a)/(p0*a+p1*b+p2*c))/length(a)
    // p1n <- 1-sum((b)/(p0*a+p1*b+p2*c))/length(a)
    map(p0, p1);

    double p2 = 1-p0-p1;
    double p0_next = 0, p1_next = 0;

    for (int i = 1; i<=commark; i++) {
      p0_next += p0*a[i]/(a[i]*p0+b[i]*p1+c[i]*p2);
      p1_next += p1*b[i]/(a[i]*p0+b[i]*p1+c[i]*p2);
    }

    nX = p0-p0_next/commark;
    nY = p1-p1_next/commark;
  }

  // The Jacobian by approximation

  void Jacob(double J[2][2], double p0, double p1)
  {
    double e = 1e-6;

    double x = p0, y = p1;
    double x_dx, x_dy, y_dx, y_dy;
    double fx, dy;

    optF(fx, fy, x, y);
    optF(x_dx, y_dx, x+e, y);
    optF(x_dy, y_dy, x, y+e);

    J[0][0] = (x_dx-fx)/e;
    J[0][1] = (x_dy-fx)/e;
    J[1][0] = (y_dx-fy)/e;
    J[1][1] = (y_dy-fy)/e;

  }

  // Analytical Jacobian. Make sure it is right!

  void Jacob2(double J[2][2], double p0, double p1)
  {
    double p0dx = 0, p0dy = 0, p1dx = 0, p1dy = 0;

    double sinp0 = sin(p0), sinp1 = sin(p1), cosp0 = cos(p0), cosp1 = cos(p1);

    for (int i = 1; i<=commark; i++) {
      double g = (a[i]-c[i])*(1+sinp0)/2+(b[i]-c[i])*(1+sinp1)*(1-sinp0)/4+c[i];

      p0dx += (0.5*a[i]*(a[i]-c[i])*cosp0-0.25*a[i]*(b[i]-c[i])*(1+sinp1)*cosp0)/(g*g);

      p0dy += (0.25*a[i]* (b[i]-c[i]) * (1-sinp0)*cosp1)/(g*g);

      p1dx += (0.5*b[i]*(a[i]-c[i])*cosp0-0.25*b[i]*(b[i]-c[i])*(1+sinp1)*cosp0)/(g*g);

      p1dy += (0.25*b[i]* (b[i]-c[i]) * (1-sinp0)*cosp1)/(g*g);
    }

    p0dx = p0dx/commark;
    p0dy = p0dy/commark;
    p1dx = p1dx/commark;
    p1dy = p1dy/commark;

    J[0][0] = p0dx;
    J[0][1] = p0dy;
    J[1][0] = p1dx;
    J[1][1] = p1dy;
  }

  void p(double J[2][2])
  {
    printf("%f %f\n%f %f\n", J[0][0], J[0][1], J[1][0], J[1][1]);
  }

  inline double max(double a, double b)
  {
    return (a>b)?a:b;
  }

  // EM algorithm using a 2 dimensional Newton-Rhapson speedup
  // in each iteration compute
  // v <- c(x,y)-qr.solve(J(x,y,e))%*%ff(x,y);

  int EM(double &X, double &Y, double eps, int miter)
  {

    X = 0.2;
    Y = 0.2;

    double x2 = X, y2 = Y, x = 99, y = 99;
    int niter = 0, countNormal = 0;

    // Initial approximation to the solution
    while(diff(x2,y2,x,y)>1e-2) { x=x2; y=y2; niter++; EM_getnext(x2,y2,x,y);  }

    // printf("Diff=%g iter=%d\n",diff(x2,y2,x,y),niter);
    // printf("******** EM ********\n");

    niter = 0;

    while (niter<miter) {
      x = x2;
      y = y2;

      double J[2][2], JI[2][2];

      Jacob(J, x, y);

      // Invert J

      double a = J[0][0], b = J[0][1], c = J[1][0], d = J[1][1];
      double det = 1/(a*d-b*c);

      if(fabs(det) < 1e-14 ) { X = -1; Y = -1; return 0; }

      JI[0][0] = d*det;
      JI[0][1] = -b*det;
      JI[1][0] = -c*det;
      JI[1][1] = a*det;

      double ffx, ffy;
      optF(ffx, ffy, x, y);

      // Newton-Rhapson iteration

      x2 = x-(JI[0][0]*ffx+JI[0][1]*ffy);
      y2 = y-(JI[1][0]*ffx+JI[1][1]*ffy);

      //printf("OPTF=%f %f   ",ffx,ffy);
      //printf("EM step: %.7f,%.7f->%.7f, %.7f | step %g %g | sum=%.7f\n",x,y,x2,y2,fabs(x2-x),fabs(y2-y), x2+y2);

      double e = max(fabs(x2-x), fabs(y2-y));

      if (e<eps) break;

      // We all pray that this worked.
      niter++;
    }

    map(x2, y2);

    if(x2 < -1e-5 || y2 < -1e-5 || (x2+y2) > 1+1e-5) { X = Y = -1; return 0; }

    if (x2+y2>1) y2 = y2+(1-x2-y2);

    if (x2<0) x2 = 0;
    if (y2<0) y2 = 0;

    X = x2;
    Y = y2;

    EM_getnext(x, y, x2, y2);
    
    //printf("EM FINISHED: (%g,%g)  steps = %d total %d newton  EPS=%g\n-----------------------\n",x2,y2,niter,niter-countNormal, fabs(x2-x) + fabs(y2-y));
    //printf("iter=%d\n",niter);

    if(diff(x2,y2,x,y) > 0.5e-6) { X = -1 ; Y = -1; return 0; }

    return 1;
  }



    #if 0
  void EM_linesearch(double &x2, double &y2, double x, double y, double eps, double miter)
  {
    double t = eps; // speedup parameter
    eps = 0.5e-6;

    double xn, yn, xq, yq;
    double d;

    int n = 1;

    while (1) {
      EM_getnext(xn, yn, x, y);
      d = diff(xn, yn, x, y);
      if (d<eps) break;

      // printf("initial x=%f y=%f D=%g dx=%g dy=%g \n",x,y,d,xn-x,yn-y);

      double x3, y3, mult = 1;

      while (1) {
        mult *= t;
        x3 = x+(xn-x)*mult;
        y3 = y+(yn-y)*mult;

        if(x3<0 || y3<0 || x3>1 || y3>1 || x3+y3 > 1) break;

        n++;
        EM_getnext(xq, yq, x3, y3);

        double d3 = diff(xq, yq, x3, y3);

        //printf("   -> %f,%f -> %f,%f mult=%f d=%g\n",x,y,x3,y3, mult, d3);

        if (d3>d) break;

        d = d3;
        //mult *=t; break;
      }

      mult /= t;

      //printf("stepped from %f,%f to %f,%f mult=%f d=%g\n",x,y,x+(xn-x)*mult, y+(yn-y)*mult ,mult,d);

      x = x+(xn-x)*mult;
      y = y+(yn-y)*mult;

    }

    x2 = x;
    y2 = y;

    //printf("Final: %f %f d=%f\n\n",x,y,d);
    printf("************* N = %d ************\n", n);
  }
  #endif

};

