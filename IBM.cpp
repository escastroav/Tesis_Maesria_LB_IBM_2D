#include "IBM.h"

IBMDisk::IBMDisk(int nDots, double r, double A0, double B0, double b, double m, double vs, double X, double Y)
{
  NDots = nDots; radius = r; A = A0; B = B0; X_center = X; Y_center = Y;
  dots_x = new double [NDots]; dots_y = new double [NDots];
  normals_x = new double [NDots];normals_y = new double [NDots];
  fx = 0; fy = 0; Vx = 0; Vy = 0; mass = m; bulk = b;
  density = mass / (M_PI*radius*radius);
  vsound = vs;
  LocateDots_Ellipse();
  LocateNormals();
}
IBMDisk::~IBMDisk(void)
{
  delete [] dots_x; delete [] dots_y;
  delete [] normals_x; delete [] normals_y;
}
double IBMDisk::Fx_J(LatticeBoltzmann & LB)
{
  double fx_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0;  
  double I_Jx = 0;
  double I_Jy = 0;
  double I_rho= 0;
  for(int k = 0; k < NDots; k++)
    {
      I_rho= LB.Interpolate('R', dots_x[k], dots_y[k]);
      I_Jx = LB.Interpolate('X', dots_x[k], dots_y[k]);
      I_Jy = LB.Interpolate('Y', dots_x[k], dots_y[k]);
      fx_tmp += - I_Jx * (I_Jx * normals_x[k] + I_Jy * normals_y[k]) * ds ;
    }
  return fx_tmp;
}
double IBMDisk::Fy_J(LatticeBoltzmann & LB)
{
  double fy_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0;  
  double I_Jx = 0;
  double I_Jy = 0;
  double I_rho= 0;
  for(int k = 0; k < NDots; k++)
    {
      I_rho = LB.Interpolate('R',dots_x[k], dots_y[k]);
      I_Jx = LB.Interpolate('X', dots_x[k], dots_y[k]);
      I_Jy = LB.Interpolate('Y', dots_x[k], dots_y[k]);
      fy_tmp += - I_Jy * (I_Jx * normals_x[k] + I_Jy * normals_y[k]) * ds ;
    }
  return fy_tmp;
}
double IBMDisk::Fx_p(LatticeBoltzmann & LB)
{
  double fx_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0; 
  double I_p = 0;
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;      
      I_p = LB.Interpolate('P', mid_dotx, mid_doty);
      fx_tmp += I_p * normals_x[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty);
  fx_tmp += - I_p * normals_x[NDots-1] * ds;
  
  return fx_tmp;
}
double IBMDisk::Fy_p(LatticeBoltzmann & LB)
{
  double fy_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0; 
  double I_p = 0;
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;      
      I_p = LB.Interpolate('P', mid_dotx, mid_doty);
      fy_tmp += I_p * normals_y[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty);
  fy_tmp += - I_p * normals_y[NDots-1] * ds;
  
  return fy_tmp;
}
double IBMDisk::Fx_p2_v2(LatticeBoltzmann & LB)
{
  double fx_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0; 
  double I_p = 0;
  double I_Jx = 0;
  double I_Jy = 0;
  double v2 = 0;
  double c2 = 0;
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;      
      I_p = LB.Interpolate('P', mid_dotx, mid_doty);
      I_Jx = LB.Interpolate('X', mid_dotx, mid_doty);
      I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty);
      v2 = I_Jx * I_Jx + I_Jy * I_Jy;
      //c2 = SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound) * SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound);
      fx_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_x[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty);
  I_Jx = LB.Interpolate('X', mid_dotx, mid_doty);
  I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty);
  v2 = I_Jx * I_Jx + I_Jy * I_Jy;
  //c2 = SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound) * SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound);
  fx_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_x[NDots-1] * ds;
  
  return fx_tmp;
}
double IBMDisk::Fy_p2_v2(LatticeBoltzmann & LB)
{
  double fy_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0;  
  double I_p = 0;
  double I_Jx = 0;
  double I_Jy = 0;
  double v2 = 0;
  double c2 = 0;
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;      
      I_p = LB.Interpolate('P', mid_dotx, mid_doty);
      I_Jx = LB.Interpolate('X', mid_dotx, mid_doty);
      I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty);
      v2 = I_Jx * I_Jx + I_Jy * I_Jy;
      //c2 = SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound) * SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound);
      fy_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_y[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty);
  I_Jx = LB.Interpolate('X', mid_dotx, mid_doty);
  I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty);
  v2 = I_Jx * I_Jx + I_Jy * I_Jy;
  //c2 = SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound) * SoundSpeed(mid_dotx,mid_doty,X_center,Y_center,radius,C,vsound);
  fy_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_y[NDots-1] * ds;

  return fy_tmp;
}
double IBMDisk::Fx_ib(LatticeBoltzmann & LB,double B, double Ux,double ds,double Rad)
{
  double fx_tmp = 0;
  for(int k=0; k<NDots;k++)
    {
      fx_tmp -= LB.Fbpx(NDots, floor(dots_x[k]), floor(dots_y[k]), dots_x, dots_y,Ux,ds,Rad) * normals_x[NDots] * ds; 
    }
  return fx_tmp;
}
double IBMDisk::Fy_ib(LatticeBoltzmann & LB,double B, double Uy,double ds,double Rad)
{
  double fy_tmp = 0;
  for(int k=0; k<NDots;k++)
    {
      fy_tmp -= LB.Fbpy(NDots, floor(dots_x[k]), floor(dots_y[k]), dots_x, dots_y,Uy,ds,Rad) * normals_y[NDots] * ds; 
    }
  return fy_tmp;
}
double IBMDisk::SoundSpeed(double x, double y, double X, double Y, double R, double c0, double c)
{
  double xp = x - X; double yp = y - Y;
  double r2 = xp*xp + yp*yp;
  double w = 1.0/3.0;
  return c0 - (c0-c)*0.5*(1-tanh(w*(r2-R*R))); 
}
void IBMDisk::MeasureForce(LatticeBoltzmann & LB, int t, int tmax, double T, double & F_min, double & F_max, double & F_add)
{
  double F_amp=Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  if(F_amp > F_max) F_max = F_amp;
  if(F_amp < F_min) F_min = F_amp;
  F_add += F_amp/T;
}
void IBMDisk::UpdatePEFRL(LatticeBoltzmann & LB, double dt)
{
  Update_X(dt,epsilon); Update_Y(dt, epsilon);

  fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,coeff1); Update_Vy(dt, coeff1);

  Update_X(dt,chi); Update_Y(dt, chi);

  fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,lambda); Update_Vy(dt, lambda);

  Update_X(dt,coeff2); Update_Y(dt, coeff2);

  fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,lambda); Update_Vy(dt, lambda);

  Update_X(dt,chi); Update_Y(dt, chi);

  fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,coeff1); Update_Vy(dt, coeff1);

  Update_X(dt,epsilon); Update_Y(dt, epsilon);
}
void IBMDisk::Print(const char * fileName, int t)
{
  ofstream MyFile(fileName);
  MyFile << t << "\t" << X_center << "\t" << Y_center << endl;
  MyFile.close();
}
