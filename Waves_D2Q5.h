#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;
#include "ComputeEpsilon.h"
//------------LB constants--------------------
const double Lambda = 500; 
const int Lx= Lambda/2 + 1;
const int Ly= Lx;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 cells
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);
const double T = Lambda / C;
const double Omega = 2*M_PI / T;
const double K = 2*M_PI / Lambda;

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
  double Speed(int ix, int iy, int X, int Y, double R, double v);
  double feq(double rho0,double Jx0,double Jy0,int i,double c);
  void Collision(int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux,double Uy, double Radius, double Ds, ComputeEpsilon & CE, double c, int t);
  void ImposeFields(int t,double X, double Y, double Radius, double c);
  void Advection(void);
  void Start(double rho0,double Jx0,double Jy0,double Fx0,double Fy0,double X, double Y, double Radius, double c);
  void Print(const char * NombreArchivo,int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux, double Uy, double Radius,double ds, ComputeEpsilon & CE);
};  

