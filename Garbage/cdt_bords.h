#ifndef CDT_BORDS_H
#define CDT_BORDS_H

typedef struct conditions_aux_bords{
  double (*f)(double, double, double);
  double (*g)(double, double, double);
  double (*h)(double, double, double);
} cdt_aux_bords;

double f1( double posx, double posy, double t);

double f2 ( double posx, double posy, double t);

double func_zero(double posx, double posy, double t);

double func_one(double posx, double posy, double t);

#endif