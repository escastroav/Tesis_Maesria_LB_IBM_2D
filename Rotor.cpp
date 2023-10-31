#include "Rotor.h"

Rotor::Rotor( double r10, double r20, double R10, double R20, double X0, double Y0, double Vx0, double Vy0, double theta0, double omega0, double rhop0, double mu0, double cp0)
{
  r1 = r10; r2 = r20;
  R1 = R10; R2 = R20;
  X_cm = X0; Y_cm = Y0;
  Vx_cm = Vx0; Vy_cm = Vy0;
  theta = theta0; omega = omega0;
  rhop = rhop0;
  mu = mu0;
  cp = cp0;
  Vol1 = (4.0*M_PI/3.0)*R1*R1*R1;
  Vol2 = (4.0*M_PI/3.0)*R2*R2*R2;
  Vol = Vol1 + Vol2;
  rhop = rhop0;
  mass1 = r2*rhop*Vol/(r1+r2); mass2 = r1*rhop*Vol/(r1+r2);
  Mass = mass1 + mass2;
  Inertia = mass1*r1*r1 + mass2*r2*r2;
  Set_dotsX(); Set_dotsY();
}
Rotor::~Rotor(void)
{
}
double Rotor::Fx_Gorkov(double c0, double rho0, double k, double P0)
{
  double rhor = rhop / rho0;
  double kappa = (c0*c0*rho0) / (cp*cp*rhop);
  double phi = (5*rhor - 2)/(2*rhor + 1) - kappa; 
  
  double Eac = P0*P0 / (4.0*rho0*c0);

  double Fx1 = Vol1 * phi * Eac * sin(2*k*x1);
  double Fx2 = Vol2 * phi * Eac * sin(2*k*x2);
  return Fx1 + Fx2;
}
double Rotor::Tz_Gorkov(double c0, double rho0, double k, double P0)
{
  double rhor = rhop / rho0;
  double kappa = (c0*c0*rho0) / (cp*cp*rhop);
  double phi = (5*rhor - 2)/(2*rhor + 1) - kappa; 
  
  double Eac = P0*P0 / (4.0*rho0*c0);

  double Fx1 = Vol1 * phi * Eac * sin(2*k*x1);
  double Fx2 = Vol2 * phi * Eac * sin(2*k*x2);

  double tau1 = - Fx1 * (y1 - Y_cm) ;
  double tau2 = - Fx2 * (y2 - Y_cm) ;
  return tau1 + tau2;
}
double Rotor::Fx_drag(double eta0)
{
  double gamma1 = 6.0 * M_PI * eta0 * R1;	
  double Fx1 = -gamma1 * (Vx_cm - omega * (y1 - Y_cm));

  double gamma2 = 6.0 * M_PI * eta0 * R2;	
  double Fx2 = -gamma2 * (Vx_cm - omega * (y2 - Y_cm));

  return Fx1 + Fx2;
}
double Rotor::Fy_drag(double eta0)
{
  double gamma1 = 6.0 * M_PI * eta0 * R1;	
  double Fy1 = -gamma1 * (Vy_cm + omega * (x1 - X_cm));

  double gamma2 = 6.0 * M_PI * eta0 * R2;	
  double Fy2 = -gamma2 * (Vy_cm + omega * (x2 - X_cm));

  return Fy1 + Fy2;
}
double Rotor::Tz_drag(double eta0)
{
  double gamma1 = 6.0 * M_PI * eta0 * R1;	
  double Fx1 = -gamma1 * (Vx_cm - omega * (y1 - Y_cm));
  double Fy1 = -gamma1 * (Vy_cm - omega * (x1 - X_cm));

  double gamma2 = 6.0 * M_PI * eta0 * R2;	
  double Fx2 = -gamma2 * (Vx_cm - omega * (y2 - Y_cm));
  double Fy2 = -gamma2 * (Vy_cm - omega * (x2 - X_cm));

  double tau1 = Fy1 * (x1 - X_cm) - Fx1 * (y1 - Y_cm);
  double tau2 = Fy2 * (x2 - X_cm) - Fx2 * (y2 - Y_cm);
  return tau1 + tau2;
}
double Rotor::Tz_Magnetic(double t, double OmegaB, double B0)
{
  return mu * B0 * sin(OmegaB*t - theta); 
}
void Rotor::UpdatePEFRL(double t, double dt, double omegaB, double c0, double rho0, double k, double P0, double eta0, double B0)
{

  Update_X(dt,epsilon); Update_Y(dt, epsilon);
  Update_Theta(dt,epsilon);
  Set_dotsX(); Set_dotsY();

  Fx_cm = Fx_Gorkov(c0, rho0, k, P0) + Fx_drag(eta0);
  Fy_cm = Fy_drag(eta0);
  Tau = Tz_Gorkov(c0, rho0, k, P0) + Tz_drag(eta0) + Tz_Magnetic(t, omegaB, B0);
  //cout << Tau << endl;

  Update_Vx(dt,coeff1); Update_Vy(dt, coeff1);
  Update_Omega(dt,coeff1);

  Update_X(dt,chi); Update_Y(dt, chi);
  Update_Theta(dt,chi);
  Set_dotsX(); Set_dotsY();

  Fx_cm = Fx_Gorkov(c0, rho0, k, P0) + Fx_drag(eta0);
  Fy_cm = Fy_drag(eta0);
  Tau = Tz_Gorkov(c0, rho0, k, P0) + Tz_drag(eta0) + Tz_Magnetic(t, omegaB, B0);

  Update_Vx(dt,lambda); Update_Vy(dt, lambda);
  Update_Omega(dt,lambda);

  Update_X(dt,coeff2); Update_Y(dt, coeff2);
  Update_Theta(dt,coeff2);
  Set_dotsX(); Set_dotsY();

  Fx_cm = Fx_Gorkov(c0, rho0, k, P0) + Fx_drag(eta0);
  Fy_cm = Fy_drag(eta0);
  Tau = Tz_Gorkov(c0, rho0, k, P0) + Tz_drag(eta0) + Tz_Magnetic(t, omegaB, B0);

  Update_Vx(dt,lambda); Update_Vy(dt, lambda);
  Update_Omega(dt,lambda);

  Update_X(dt,chi); Update_Y(dt, chi);
  Update_Theta(dt,chi);
  Set_dotsX(); Set_dotsY();

  Fx_cm = Fx_Gorkov(c0, rho0, k, P0) + Fx_drag(eta0);
  Fy_cm = Fy_drag(eta0);
  Tau = Tz_Gorkov(c0, rho0, k, P0) + Tz_drag(eta0) + Tz_Magnetic(t, omegaB, B0);

  Update_Vx(dt,lambda); Update_Vy(dt, lambda);
  Update_Omega(dt,lambda);

  Update_X(dt,epsilon); Update_Y(dt, epsilon);
  Update_Theta(dt,epsilon);
  Set_dotsX(); Set_dotsY();


}
