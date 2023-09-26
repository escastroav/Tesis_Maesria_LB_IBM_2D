#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

//----------MD constants-----------------------
const double epsilon = 0.1786178958448091;
const double lambda = -0.2123418310626054;
const double chi = -0.06626458266981849;

const double coeff1 = (1 - 2*lambda)*0.5;
const double coeff2 = (1-2*(chi+epsilon));
//--------------------- Molecular dynamics----------------
class Rotor{
private:
  double radius;
  double r1, r2;
  double Vol, Vol1, Vol2;
  double rhop;
  double cp;
  double X_cm; double Y_cm;
  double x1, x2;
  double y1, y2;
  double Vx_cm; double Vy_cm;
  double theta, omega;
  double mu;
  // Molecular dynamics
  double mass1, mass2;
  double Fx_cm, Fy_cm;
  double Tau;
  double Mass, Inertia;
  void Set_Rs(){r1 = radius; r2 = mass2 * radius / mass1;};
  void Update_Vx(double dt, double coeff){Vx_cm+=Fx_cm*(dt*coeff/Mass);};
  void Update_Vy(double dt, double coeff){Vy_cm+=Fy_cm*(dt*coeff/Mass);};
  void Update_Omega(double dt, double coeff){omega+=Tau*(dt*coeff/Inertia);};
  void Update_X(double dt, double coeff){X_cm += Vx*(dt*coeff); x1 = X_cm - r1 * cos(theta); x2 = X_cm + r2 * cos(theta);};
  void Update_Y(double dt, double coeff){Y_cm += Vy*(dt*coeff); y1 = Y_cm - r1 * sin(theta); y2 = Y_cm + r2 * sin(theta);};
  void Update_Theta(double dt, double coeff){theta+=omega*(dt*coeff);};
public:
  Rotor(double r0, double X0, double Y0, double theta0, double omega0, double m1, double m2 ,double mu0);
  ~Rotor(void);
  //Getters
  double GetMass(){return Mass;}; 
  double GetMass1(){return mass1;}; 
  double GetMass2(){return mass2;};
  double GetInertia(){return Inertia;};
  double GetRadius(){return radius;};
  double GetGamma1(){return gamma1;};
  double GetGamma2(){return gamma2;};
  double GetX(){return X_cm;};
  double GetY(){return Y_cm;};
  double GetX1(){return x1;};double GetX2(){return x2;};
  double GetY1(){return y1;};double GetY2(){return y2;};
  double GetTheta(){return theta;};
  double GetOmega(){return omega;};
  double GetVx(){return Vx_cm;};
  double GetVy(){return Vy_cm;};
  double GetFx(){return Fx_cm;};
  double GetFy(){return Fy_cm;};
  double GetTau(){return Tau;};
  //Forces applied to center of mass
  double Fx_Gorkov(double c0, double rho0, double k, double P0);
  double Fx_drag(double eta0);
  double Fy_drag(double eta0);
  //Torques applied around center of mass
  double Tz_Gorkov(double c0, double rho0, double k, double P0);
  double Tz_Magnetic(double t, double omegaB, double B0);
  double Tz_drag(double eta0);
  //Integrator
  void UpdatePEFRL(double t, double dt, double omegaB, double c0, double rho0, double k, double P0, double eta0, double B0);
};
