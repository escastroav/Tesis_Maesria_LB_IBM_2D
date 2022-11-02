#include "IBM.h"

IBMDisk::IBMDisk(int nDots, double rad, double X, double Y)
{
  NDots = nDots; radius = rad; X_center = X; Y_center = Y;
  dots_x = new double [NDots]; dots_y = new double [NDots];
  normals_x = new double [NDots];normals_y = new double [NDots];
  fx = 0; fy = 0; Vx = 0; Vy = 0; mass = 10.0; bulk = 0.5;
  LocateDots();
}
IBMDisk::~IBMDisk(void)
{
  delete [] dots_x; delete [] dots_y;
  delete [] normals_x; delete [] normals_y;
}
void IBMDisk::Fx(LatticeBoltzmann & LB)
{
  double fx_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0;  
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;      
      fx_tmp += LB.Interpolate('P', mid_dotx, mid_doty) * normals_x[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  fx_tmp += LB.Interpolate('P', mid_dotx, mid_doty) * normals_x[NDots-1] * ds;
  
  fx = fx_tmp;
}
void IBMDisk::Fy(LatticeBoltzmann & LB)
{
  double fy_tmp = 0;
  double mid_dotx = 0;
  double mid_doty = 0;  
  for(int k = 0; k < NDots-1; k++)
    {
      mid_dotx = (dots_x[k+1] + dots_x[k]) * 0.5;
      mid_doty = (dots_y[k+1] + dots_y[k]) * 0.5;      
      fy_tmp += LB.Interpolate('P', mid_dotx, mid_doty) * normals_y[k] * ds;
    }
  mid_dotx = (dots_x[0] + dots_x[NDots-1]) * 0.5;
  mid_doty = (dots_y[0] + dots_y[NDots-1]) * 0.5;      
  fy_tmp += LB.Interpolate('P', mid_dotx, mid_doty) * normals_y[NDots-1] * ds;

  fy = fy_tmp;
}
void IBMDisk::UpdatePEFRL(LatticeBoltzmann & LB, double dt)
{
  Update_X(dt,epsilon); Update_Y(dt, epsilon);

  Fx(LB); Fy(LB);
  Update_Vx(dt,coeff1); Update_Vy(dt, coeff1);

  Update_X(dt,chi); Update_Y(dt, chi);

  Fx(LB); Fy(LB);
  Update_Vx(dt,lambda); Update_Vy(dt, lambda);

  Update_X(dt,coeff2); Update_Y(dt, coeff2);

  Fx(LB); Fy(LB);
  Update_Vx(dt,lambda); Update_Vy(dt, lambda);

  Update_X(dt,chi); Update_Y(dt, chi);

  Fx(LB); Fy(LB);
  Update_Vx(dt,coeff1); Update_Vy(dt, coeff1);

  Update_X(dt,epsilon); Update_Y(dt, epsilon);
}
void IBMDisk::Print(const char * fileName, int t)
{
  ofstream MyFile(fileName);
  MyFile << t << "\t" << X_center << "\t" << Y_center << endl;
  MyFile.close();
}
