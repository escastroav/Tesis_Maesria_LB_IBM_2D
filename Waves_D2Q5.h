#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;
#include "ComputeEpsilon.h"
//------------LB constants--------------------
const double Lambda = 500; 
const int Lx= 250;
const int Ly= 64; 
const int Lz= 64;

const int Q=7;
const double W0=1.0/4;

const double C=0.25;
const double C2=C*C;
const double AUX0=1-4*C2*(1-W0);
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
  int Vx[Q],Vy[Q],Vz[Q];  //Velocity vectors
  double *f, *fnew; //Distribution Functions
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int iz,int i){return ((iz*Ly+iy)*Lx+ix)*Q+i;};
  double rho(int ix,int iy,int iz,bool UseNew);
  double Jx(int ix,int iy,int iz,bool UseNew);
  double Jy(int ix,int iy,int iz,bool UseNew);
  double Jz(int ix,int iy,int iz,bool UseNew);
  double Pi(int ix,int iy,int iz,bool UseNew);
  //Interpolation
  double Interpolate(char field, double x, double y, double z);
  //Eq distribution
  double Speed(int ix, int iy, int iz, int X, int Y, int Z, double R, double v);
  double feq(double rho0,double Jx0,double Jy0,double Jz0,int i,double c);
  void Collision(double X, double Y, double Z, double Radius, double c);
  void ImposeFields(int t,double X,double Y,double Z, double Radius, double c, double Po);
  void Advection(void);
  void Start(double rho0,double Jx0,double Jy0,double Jz0,double X, double Y, double Z, double Radius, double c);
  void PrintBoundary(const char * NameFile, int Ndots, double * dotsx, double * dotsy, double X, double Y);
  void Print(const char * NombreArchivo, double X, double Y, double Z, double Radius,double v);
};  

