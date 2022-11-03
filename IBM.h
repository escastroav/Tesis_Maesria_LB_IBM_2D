#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

#include "Waves_D2Q5.h"
//----------MD constants-----------------------
const double epsilon = 0.1786178958448091;
const double lambda = -0.2123418310626054;
const double chi = -0.06626458266981849;

const double coeff1 = (1 - 2*lambda)*0.5;
const double coeff2 = (1-2*(chi+epsilon));
//--------------------- ImmerseBoundary & Molecular dynamics----------------
class IBMDisk{
private:
  int NDots;
  double radius;
  double ds;
  double X_center; double Y_center;
  double *dots_x;
  double *dots_y;
  double *normals_x;
  double *normals_y;
  // Molecular dynamics
  double mass; double bulk; double density;
  double fx; double fy;
  double Vx; double Vy;
  void Update_Vx(double dt, double coeff){Vx+=fx*(dt*coeff/mass);};
  void Update_Vy(double dt, double coeff){Vy+=fy*(dt*coeff/mass);};
  void Update_X(double dt, double coeff)
  {
    bool isOut = false;
    for(int k=0; k<NDots; k++){
      isOut = !(dots_x[k] > 0 && dots_x[k] < Lx);
      if(isOut){
	//cout << "Ball is out!" << endl;
	return;}
      else
	dots_x[k] += Vx*(dt*coeff);
    }
    if(!isOut)
      X_center += Vx*(dt*coeff);
  };
  void Update_Y(double dt, double coeff)
  {
    bool isOut = false;
    for(int k=0; k<NDots; k++){
      isOut = !(dots_y[k] > 0 && dots_y[k] < Ly);
      if(isOut){
        //cout << "Ball is out!" << endl;
	return;}
      else
	dots_y[k] += Vy*(dt*coeff);
    }
    if(!isOut)
      Y_center += Vy*(dt*coeff);
  };
  void LocateDots(void)
  {
    int k = 0; double dtheta = 2.0*M_PI / NDots;
    ds = radius * dtheta;
    for(k = 0; k < NDots; k++)
      {
	dots_x[k] = X_center + radius * cos(dtheta * k);
	normals_x[k] = (-1.0) * cos(dtheta * (k + 0.5));
	dots_y[k] = Y_center + radius * sin(dtheta * k);
	normals_y[k] = (-1.0) * sin(dtheta * (k + 0.5));
      }
  }
public:
  IBMDisk(int nDots, double r, double b, double m, double X, double Y);
  ~IBMDisk(void);
  int GetNdots(){return NDots;};
  double GetMass(){return mass;}; double GetBulk(){return bulk;};
  double GetDensity(){return density;};
  double* GetDotsX(){return dots_x;};
  double* GetDotsY(){return dots_y;};
  double GetX(){return X_center;};
  double GetY(){return Y_center;};
  double GetVx(){return Vx;};
  double GetVy(){return Vy;};
  double GetFx(){return fx;};
  double GetFy(){return fy;};
  //Force
  void Fx(LatticeBoltzmann & LB);
  void Fy(LatticeBoltzmann & LB);
  //Integrator
  void UpdatePEFRL(LatticeBoltzmann & LB, double dt);
  //PrintDots
  void Print(const char * fileName, int t);
};
