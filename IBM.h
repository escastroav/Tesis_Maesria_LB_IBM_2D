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
  double gamma;
  double sound_speed;
  double A, B;
  double ds;
  double X_center; double Y_center;
  double phi;
  double *dots_x;
  double *dots_y;
  double *normals_x;
  double *normals_y;
  // Molecular dynamics
  double mass; double bulk; double density; double vsound;
  double fx; double fy;
  double fx_J; double fy_J;
  double Vx; double Vy;
  void Update_Vx(double dt, double coeff){Vx+=fx*(dt*coeff/mass);};
  void Update_Vy(double dt, double coeff){Vy+=fy*(dt*coeff/mass);};
  void Update_X(double dt, double coeff)
  {
    bool isOut = false;
    for(int k=0; k<NDots; k++){
      isOut = !(dots_x[k] > 0 && dots_x[k] < Lx);
      if(isOut){
	cout << "Ball is out!" << endl;
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
  void CorrectX(double dx)
  {
    bool isOut = false;
    for(int k=0; k<NDots; k++){
      isOut = !(dots_x[k] > 0 && dots_x[k] < Lx);
      if(isOut){
	cout << "Ball is out!" << endl;
	return;}
      else
	dots_x[k] += dx;
    }
    if(!isOut)
      X_center += dx;
  };
  void CorrectY(double dy)
  {
    bool isOut = false;
    for(int k=0; k<NDots; k++){
      isOut = !(dots_y[k] > 0 && dots_y[k] < Ly);
      if(isOut){
        //cout << "Ball is out!" << endl;
	return;}
      else
	dots_y[k] += dy;
    }
    if(!isOut)
      Y_center += dy;
  };
  void Rotate(void)
  {
    double rel_x, rel_y, rot_x, rot_y;
    for(int k=0; k<NDots; k++){
      rel_x = dots_x[k] - X_center;
      rel_y = dots_y[k] - Y_center;
      rot_x = rel_x * cos(phi) - rel_y * sin(phi);
      rot_y = rel_y * cos(phi) + rel_x * sin(phi);
      dots_x[k] = X_center + rot_x;
      dots_y[k] = Y_center + rot_y;
      //cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;
    } 
  };
  void LocateDots_Circle(void)
  {
    int k = 0; double dtheta = 2.0*M_PI / NDots;
    ds = radius * dtheta;
    for(k = 0; k < NDots; k++)
      {
	dots_x[k] = X_center + radius * cos(dtheta * k);
	//normals_x[k] = cos(dtheta * ((double)k + 0.5));
	dots_y[k] = Y_center + radius * sin(dtheta * k);
	//normals_y[k] = sin(dtheta * ((double)k + 0.5));
	//cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;
	//cout << k << " " << normals_x[k] << " " << normals_y[k] << endl;
      }    
  }
  void LocateDots_Capsule(void)
  {
    int k = 0; double dtheta = 4 * M_PI / NDots;
    ds = radius * dtheta;
    for(k = 0; k < NDots/4; k++)
      {
	dots_x[k] = X_center + radius + (2 * radius / M_PI)  * cos(dtheta * k - M_PI * 0.5);
	//normals_x[k] = cos(dtheta * ((double)k + 0.5));
	dots_y[k] = Y_center + (2 * radius / M_PI) * sin(dtheta * k - M_PI * 0.5);
	//normals_y[k] = sin(dtheta * ((double)k + 0.5));
	//cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;
	//cout << k << " " << normals_x[k] << " " << normals_y[k] << endl;
      }    
    for(k = NDots/4; k < NDots/2; k++)
      {
	dots_x[k] = X_center + radius - k + NDots * 0.25;
	//normals_x[k] = cos(dtheta * ((double)k + 0.5));
	dots_y[k] = Y_center + (2 * radius / M_PI);
	//normals_y[k] = sin(dtheta * ((double)k + 0.5));
	//cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;
	//cout << k << " " << normals_x[k] << " " << normals_y[k] << endl;
      }    
    for(k = NDots/2; k < 3*NDots/4; k++)
      {
	dots_x[k] = X_center - radius + (2 * radius / M_PI) * cos(dtheta * (k - NDots/4) - M_PI * 0.5);
	//normals_x[k] = cos(dtheta * ((double)k + 0.5));
	dots_y[k] = Y_center + (2 * radius / M_PI ) * sin(dtheta * (k - NDots/4) - M_PI * 0.5);
	//normals_y[k] = sin(dtheta * ((double)k + 0.5));
	//cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;
	//cout << k << " " << normals_x[k] << " " << normals_y[k] << endl;
      }    
    for(k = 3*NDots/4; k < NDots; k++)
      {
	dots_x[k] = X_center - radius + k - NDots * 0.75;
	//normals_x[k] = cos(dtheta * ((double)k + 0.5));
	dots_y[k] = Y_center - (2 * radius / M_PI);
	//normals_y[k] = sin(dtheta * ((double)k + 0.5));
	//cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;
	//cout << k << " " << normals_x[k] << " " << normals_y[k] << endl;
      } 
    /**for(k = 0; k < NDots; k++)
    {   
	cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;
    }**/
  }
  void LocateDots_Ellipse(void)
  {
    int k = 0; double dtheta = 2.0*M_PI / NDots;
    ds = radius * dtheta;
    for(k = 0; k < NDots; k++)
      {
	dots_x[k] = X_center + A * cos(dtheta * k);// * cos(phi) - B * sin(dtheta * k) * sin(phi);
	dots_y[k] = Y_center + B * sin(dtheta * k);// * cos(phi) + A * cos(dtheta * k) * sin(phi);
	//cout << k << " " << dots_x[k] << " " << dots_y[k] << endl;

      }    
  }
  void LocateNormals(void)
  {
    int k = 0; double dtheta = 2.0*M_PI / NDots;
    double diff_x, diff_y, diff;
    ds = radius * dtheta;
    for(k = 0; k < NDots-1; k++)
      {
	diff_x = dots_x[k+1] - dots_x[k];
	diff_y = dots_y[k+1] - dots_y[k];
	diff = sqrt(diff_x*diff_x + diff_y*diff_y);
	normals_x[k] = diff_y / diff;
	normals_y[k] = -diff_x / diff;
	//cout << k << " " << normals_x[k] << " " << normals_y[k] << endl;
	//cout << k << " " << (dots_x[k]+dots_x[k+1])*0.5 << " " << normals_x[k] << " " << (dots_y[k]+dots_y[k+1])*0.5<< " " << normals_y[k] << endl;
      }
    diff_x = dots_x[0] - dots_x[NDots - 1];
    diff_y = dots_y[0] - dots_y[NDots - 1];
    diff = sqrt(diff_x*diff_x + diff_y*diff_y);
    normals_x[NDots - 1] = diff_y / diff;
    normals_y[NDots - 1] = -diff_x / diff;
	//cout << k << " " << (dots_x[NDots-1]+dots_x[0])*0.5 << " " << normals_x[k] << " " << (dots_y[NDots-1]+dots_y[0])*0.5<< " " << normals_y[k] << endl;
    //cout << k << " " << normals_x[NDots-1] << " " << normals_y[NDots - 1] << endl;
  }
public:
  IBMDisk(int nDots, double r, double g, double Phi0, double A0, double B0, double b, double m, double vs, double X, double Y);
  ~IBMDisk(void);
  int GetNdots(){return NDots;}; double GetDs(){return ds;};
  double GetMass(){return mass;}; double GetBulk(){return bulk;};
  double GetSoundSpeed(){return vsound;};
  double GetRadius(){return radius;};
  double GetDensity(){return density;};
  double* GetDotsX(){return dots_x;};
  double* GetDotsY(){return dots_y;};
  double GetX(){return X_center;};
  double GetY(){return Y_center;};
  double GetPhi(){return phi;};
  double GetVx(){return Vx;};
  double GetVy(){return Vy;};
  double GetFx(){return fx;};
  double GetFx_J(){return fx_J;};
  double GetFy_J(){return fy_J;};
  double GetFy(){return fy;};
  double GetA(){return A;};
  double GetB(){return B;};
  //Forces
  double Fx(LatticeBoltzmann & LB);
  double Fx_p(LatticeBoltzmann & LB);
  double Fx_p2_v2(LatticeBoltzmann & LB);
  double Fx_J(LatticeBoltzmann & LB);
  double Fy(LatticeBoltzmann & LB);
  double Fy_p(LatticeBoltzmann & LB);
  double Fy_p2_v2(LatticeBoltzmann & LB);
  double Fy_J(LatticeBoltzmann & LB);
  double Fx_drag();
  double Fy_drag();
  double Fx_spring(double diff, double DX, double DY);
  double Fy_spring(double diff, double DX, double DY);
  double Fx_magnetic(double t, double OmegaB, double DX, double DY);
  double Fy_magnetic(double t, double OmegaB, double DX, double DY);
  //Torque
  double Tz(LatticeBoltzmann & LB,double d, double Phi0);
  double Tz_p(LatticeBoltzmann & LB);
  double Tz_p2_v2(LatticeBoltzmann & LB);
  double Tz_J(LatticeBoltzmann & LB);
  //Integrator
  void UpdatePEFRL(double t, double omega, double ARF_x, double ARF_y, double dt,IBMDisk & Disk,double sign,double D);//(LatticeBoltzmann & LB, double dt);
  //PrintDots
  void Print(const char * fileName, int t);
};
