#include "cdt_bords.h"
#include "tools.h"
// equivalent à module fonctions.
#include <math.h>

// Résolution 3.
double f1( double posx, double posy, double t)
{
  UNUSED(t);
  double function;
  function = sin(posx) + cos(posy);
  return(function);
}

// Résolution 2.
double f2 ( double posx, double posy, double t)
{
  UNUSED(t);
  return 2*(posy - posy*posy + posx - posx*posx);
}


// par des 0.
double func_zero(double posx, double posy, double t)
{
  UNUSED(t);
  UNUSED(posx);
  UNUSED(posy); 
  return 0;
}


// par des 1.
double func_one(double posx, double posy, double t)
{
  UNUSED(t);
  UNUSED(posx);
  UNUSED(posy);
  return 1;
}
