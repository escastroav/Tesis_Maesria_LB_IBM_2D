#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;
#include "ComputeEpsilon.h"
//------------LB constants--------------------
const int Lx=128;
const int Ly=128;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//--------------------- class LatticeBoltzmann ----------------------------------
class LatticeBoltzmann{
private:
  double w[Q];      //Weights 
  int Vx[Q],Vy[Q];  //Velocity vectors
  double *f, *fnew; //Distribution Functions
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  double rho(int ix,int iy,bool UseNew);
  double Fbpx(int Ndots, int ix, int iy, double * dotsx, double * dotsy, double bulk, double Ux, double Rad, ComputeEpsilon & CE);
  double Fbpy(int Ndots, int ix, int iy, double * dotsx, double * dotsy, double bulk, double Uy, double Rad, ComputeEpsilon & CE);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  //Interpolation
  double Interpolate(char field, double x, double y);
  //Eq distribution
  double feq(double rho0,double Jx0,double Jy0,int i);
  void Collision(int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux,double Uy, double Radius, double Ds, ComputeEpsilon & CE);
  void ImposeFields(int t);
  void Advection(void);
  void Start(double rho0,double Jx0,double Jy0,double Fx0,double Fy0);
  void Print(const char * NombreArchivo,int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux, double Uy, double Radius,double ds, ComputeEpsilon & CE);
};  

