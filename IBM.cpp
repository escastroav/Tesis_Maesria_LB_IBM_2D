#include "IBM.h"

IBMDisk::IBMDisk(int nDots, double r, double g, double phi0, double A0, double B0, double b, double m, double vs, double X, double Y)
{
  NDots = nDots; radius = r; gamma = g; phi = phi0; A = A0; B = B0; X_center = X; Y_center = Y;
  dots_x = new double [NDots]; dots_y = new double [NDots];
  normals_x = new double [NDots];normals_y = new double [NDots];
  fx = 0; fy = 0; Vx = 0; Vy = 0; mass = m; bulk = b;
  density = mass / (M_PI*radius*radius);
  vsound = vs;
  LocateDots_Circle();
  //LocateDots_Ellipse();
  //Rotate();
  //LocateNormals();
}
IBMDisk::~IBMDisk(void)
{
  delete [] dots_x; delete [] dots_y;
  delete [] normals_x; delete [] normals_y;
}
double IBMDisk::Fx(LatticeBoltzmann & LB)
{
  double fx_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0;  
  double I_Jx = 0;
  double I_Jy = 0;
  double v2 = 0;
  double I_p= 0;
  double fx_J = 0;
  double fx_v2_p2= 0;
  double fx_p = 0;
  for(int k = 0; k < NDots; k++)
    {
      I_p  = LB.Interpolate('P', dots_x[k], dots_y[k], dots_z[k]);
      I_Jx = LB.Interpolate('X', dots_x[k], dots_y[k], dots_z[k]);
      I_Jy = LB.Interpolate('Y', dots_x[k], dots_y[k], dots_z[k]);
      
      v2 = I_Jx * I_Jx + I_Jy * I_Jy;
      fx_p = - I_p * normals_x[k] * ds;
      fx_J = - I_Jx * (I_Jx * normals_x[k] + I_Jy * normals_y[k]) * ds;
      fx_v2_p2 = 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_x[k] * ds;
      fx_tmp += fx_v2_p2 + fx_J;
    }
  return fx_tmp;
}
double IBMDisk::Fy(LatticeBoltzmann & LB)
{
  double fy_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0;  
  double I_Jx = 0;
  double I_Jy = 0;
  double v2 = 0;
  double I_p= 0;
  double fy_J = 0;
  double fy_v2_p2= 0;
  double fy_p = 0;
  for(int k = 0; k < NDots; k++)
    {
      I_p  = LB.Interpolate('P', dots_x[k], dots_y[k], dots_z[k]);
      I_Jx = LB.Interpolate('X', dots_x[k], dots_y[k], dots_z[k]);
      I_Jy = LB.Interpolate('Y', dots_x[k], dots_y[k], dots_z[k]);
      
      v2 = I_Jx * I_Jx + I_Jy * I_Jy;
      fy_p = - I_p * normals_x[k] * ds;
      fy_J = - I_Jx * (I_Jx * normals_x[k] + I_Jy * normals_y[k]) * ds ;
      fy_v2_p2 = 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_x[k] * ds;
      fy_tmp += fy_v2_p2 + fy_J;
    }
  return fy_tmp;
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
      I_rho= LB.Interpolate('R', dots_x[k], dots_y[k], dots_z[k]);
      I_Jx = LB.Interpolate('X', dots_x[k], dots_y[k], dots_z[k]);
      I_Jy = LB.Interpolate('Y', dots_x[k], dots_y[k], dots_z[k]);
      fx_tmp += - I_Jx * (I_Jx * normals_x[k] + I_Jy * normals_y[k]) * ds;
    }
  return fx_tmp;
}
double IBMDisk::Tz_J(LatticeBoltzmann & LB)
{
  double tz_tmp = 0;
  double rot_dotx = 0;
  double rot_doty = 0;
  double I_Jx = 0;
  double I_Jy = 0;
  double I_rho= 0;
  for(int k = 0; k < NDots; k++)
    {
      rot_dotx = dots_x[k] - X_center;
      rot_doty = dots_y[k] - Y_center;
      I_rho= LB.Interpolate('R', dots_x[k], dots_y[k], dots_z[k]);
      I_Jx = LB.Interpolate('X', dots_x[k], dots_y[k], dots_z[k]);
      I_Jy = LB.Interpolate('Y', dots_x[k], dots_y[k], dots_z[k]);
      tz_tmp += - ( rot_dotx * I_Jy - rot_doty * I_Jx ) * (I_Jx * normals_x[k] + I_Jy * normals_y[k]) * ds; 
      //cout <<  ( rot_dotx * I_Jy - rot_doty * I_Jx ) << endl;
    }
  return tz_tmp;
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
      I_rho = LB.Interpolate('R',dots_x[k], dots_y[k], dots_z[k]);
      I_Jx = LB.Interpolate('X', dots_x[k], dots_y[k], dots_z[k]);
      I_Jy = LB.Interpolate('Y', dots_x[k], dots_y[k], dots_z[k]);
      fy_tmp += - I_Jy * (I_Jx * normals_x[k] + I_Jy * normals_y[k]) * ds ;
    }
  return fy_tmp;
}
double IBMDisk::Fx_p(LatticeBoltzmann & LB)
{
  double fx_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0; 
  double mid_dotz = 0;
  double I_p = 0;
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;      
      I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_dotz);
      fx_tmp += - I_p * normals_x[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_dotz);
  fx_tmp += - I_p * normals_x[NDots-1] * ds;
  
  return fx_tmp;
}
double IBMDisk::Tz_p(LatticeBoltzmann & LB)
{
  double tz_tmp = 0;
  double mid_dotx = 0;
  double rot_dotx = 0;
  double mid_doty = 0; 
  double rot_doty = 0; 
  double I_p = 0;
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;
      rot_dotx = dots_x[k] - X_center;      
      rot_doty = dots_y[k] - Y_center;      
      I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
      tz_tmp += - I_p * ( rot_dotx * normals_y[k] - rot_doty * normals_x[k] )* ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  rot_dotx = mid_dotx - X_center;      
  rot_doty = mid_doty - Y_center;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
  tz_tmp += - I_p * ( rot_dotx * normals_y[NDots-1] - rot_doty * normals_x[NDots-1] ) * ds;
  
  return tz_tmp;
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
      I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
      fy_tmp += - I_p * normals_y[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
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
      I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
      I_Jx = LB.Interpolate('X', mid_dotx, mid_doty, mid_doty);
      I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty, mid_doty);
      v2 = I_Jx * I_Jx + I_Jy * I_Jy;
      fx_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_x[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
  I_Jx = LB.Interpolate('X', mid_dotx, mid_doty, mid_doty);
  I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty, mid_doty);
  v2 = I_Jx * I_Jx + I_Jy * I_Jy;
  fx_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_x[NDots-1] * ds;
  
  return fx_tmp;
}
double IBMDisk::Tz_p2_v2(LatticeBoltzmann & LB)
{
  double tz_tmp = 0;
  double mid_dotx = 0;
  double rot_dotx = 0;
  double mid_doty = 0;
  double rot_doty = 0; 
  double I_p = 0;
  double I_Jx = 0;
  double I_Jy = 0;
  double v2 = 0;
  double c2 = 0;
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5; 
      rot_dotx = dots_x[k] - X_center;
      rot_doty = dots_y[k] - Y_center;     
      I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
      I_Jx = LB.Interpolate('X', mid_dotx, mid_doty, mid_doty);
      I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty, mid_doty);
      v2 = I_Jx * I_Jx + I_Jy * I_Jy;
      tz_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * ( rot_dotx * normals_y[k] - rot_doty * normals_x[k] ) * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  rot_dotx = mid_dotx - X_center;
  rot_doty = mid_doty - Y_center;     
  I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
  I_Jx = LB.Interpolate('X', mid_dotx, mid_doty, mid_doty);
  I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty, mid_doty);
  v2 = I_Jx * I_Jx + I_Jy * I_Jy;
  tz_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * ( rot_dotx * normals_y[NDots-1] - rot_doty * normals_x[NDots-1] ) * ds;
  
  return tz_tmp;
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
      I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
      I_Jx = LB.Interpolate('X', mid_dotx, mid_doty, mid_doty);
      I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty, mid_doty);
      v2 = I_Jx * I_Jx + I_Jy * I_Jy;
      fy_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_y[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  I_p = LB.Interpolate('P', mid_dotx, mid_doty, mid_doty);
  I_Jx = LB.Interpolate('X', mid_dotx, mid_doty, mid_doty);
  I_Jy = LB.Interpolate('Y', mid_dotx, mid_doty, mid_doty);
  v2 = I_Jx * I_Jx + I_Jy * I_Jy;
  fy_tmp += 0.5 * ( v2 - (I_p * I_p / C2) ) * normals_y[NDots-1] * ds;

  return fy_tmp;
}
double IBMDisk::Fx_drag()
{
  return -1.0 * gamma * radius * Vx;
}
double IBMDisk::Fy_drag()
{
  return -1.0 * gamma * radius * Vy;
}
void IBMDisk::UpdatePEFRL(double ARF_x, double ARF_y, double dt)//(LatticeBoltzmann & LB, double dt)
{
  Update_X(dt,epsilon); Update_Y(dt, epsilon);

  fx = ARF_x + Fx_drag();//fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = ARF_y + Fy_drag();//fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,coeff1); Update_Vy(dt, coeff1);

  Update_X(dt,chi); Update_Y(dt, chi);

  fx = ARF_x + Fx_drag();//fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = ARF_y + Fy_drag();//fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,lambda); Update_Vy(dt, lambda);

  Update_X(dt,coeff2); Update_Y(dt, coeff2);

  fx = ARF_x + Fx_drag();//fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = ARF_y + Fy_drag();//fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,lambda); Update_Vy(dt, lambda);

  Update_X(dt,chi); Update_Y(dt, chi);

  fx = ARF_x + Fx_drag();//fx = Fx_p(LB)+Fx_p2_v2(LB)+Fx_J(LB);
  fy = ARF_y + Fy_drag();//fy = Fy_p(LB)+Fy_p2_v2(LB)+Fy_J(LB);
  Update_Vx(dt,coeff1); Update_Vy(dt, coeff1);

  Update_X(dt,epsilon); Update_Y(dt, epsilon);
}
void IBMDisk::Print(const char * fileName, int t)
{
  ofstream MyFile(fileName);
  MyFile << t << "\t" << X_center << "\t" << Y_center << endl;
  MyFile.close();
}
