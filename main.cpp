#include "Waves_D2Q5.h"
#include "IBM.h"
#include "ComputeEpsilon.h"
#include <string>
//---------------------global Interpolation functions---------------------------
double Kernel(double r)
{
  // r = x - x1 / (x0 - x1)
  if(abs(r) < 0.5)
    return (1.0+sqrt(1.0-3.0*r*r))/3.0;
  else if(abs(r) >= 0.5 && abs(r) <= 1.5)
    return (5.0-3.0*abs(r)-sqrt(1.0 - 3.0*(1.0-abs(r))*(1.0-abs(r)) ))/6.0;
  else
    return 0;
}
void Neighbours(double* nb_x, double* nb_y, double x, double y)
{
  /*
  7	8	9	10
    
  6	1	2	11
    
  5	0	3	12
    
  4	15	14	13
  */
  //first order:
  nb_x[0] = floor(x);	nb_y[0] = floor(y);
  nb_x[1] = nb_x[0];	nb_y[1] = nb_y[0]+1;
  nb_x[2] = nb_x[0]+1;	nb_y[2] = nb_y[0]+1;
  nb_x[3] = nb_x[0]+1;	nb_y[3] = nb_y[0];
  //second order:
  nb_x[4] = floor(x)-1;	nb_y[4] = floor(y)-1;
  //left sied:
  nb_x[5] = nb_x[4];	nb_y[5] = nb_y[4]+1;
  nb_x[6] = nb_x[4];	nb_y[6] = nb_y[4]+2;
  nb_x[7] = nb_x[4];	nb_y[7] = nb_y[4]+3;
  //up side:
  nb_x[8] = nb_x[4]+1;	nb_y[8] = nb_y[4]+3;
  nb_x[9] = nb_x[4]+2;	nb_y[9] = nb_y[4]+3;
  nb_x[10] = nb_x[4]+3;	nb_y[10] = nb_y[4]+3;
  //right side:
  nb_x[11] = nb_x[4]+3;	nb_y[11] = nb_y[4]+2;
  nb_x[12] = nb_x[4]+3;	nb_y[12] = nb_y[4]+1;
  nb_x[13] = nb_x[4]+3;	nb_y[13] = nb_y[4];
  //down side: 
  nb_x[14] = nb_x[4]+2;	nb_y[14] = nb_y[4];
  nb_x[15] = nb_x[4]+1;	nb_y[15] = nb_y[4];
}


//--------------- Global Functions ------------

int main(int argc, char * argv[]){
  //if(!(argc-1)){cout << "No timesteps!" << endl; return 0;}	
  LatticeBoltzmann Waves;
  int Nodes =120;
  double R = Nodes / (2*M_PI);
   double mass = 100;
  double vs = 0.3;//C=0.5
  double density = mass / (M_PI*R*R);
  double B = vs*vs*density;
  double F_max = 0;
  double F_amp = 0;
  int X0 = Lx/2;
  int Y0 = Ly/2;
  ComputeEpsilon EpSolver(Nodes);
  IBMDisk Disk(Nodes, R, B, mass, X0, Y0);
  int t,tmax=8*T;//atoi(argv[1]);
  double rho0=0,Jx0=0,Jy0=0;
  string outName = "Waves3D_st=";
  string extension = ".dat";
  string frame = "";
  //double p_max[Lx];
  //double p_min[Lx]; 
  //double pix = 0;
  //double T = 97 / C;
  int iy =  Ly/2; int ix = 0;
  //Rotores magnéticos moviéndose en paredes acústicas: Un estudio computacional
  EpSolver.BuildMatrix(Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetDs()); EpSolver.SolveA();
  //Start
  Waves.Start(rho0,Jx0,Jy0,0,0,X0,Y0,R,vs);
  /**for(int ix=0;ix<Lx;ix++)
  {
	  p_max[ix] = p_min[ix] = 0;
  }**/
  //Run
  for(t=0;t<tmax;t++){
    Waves.Collision(Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R, Disk.GetDs(), EpSolver, vs, t);
    //Waves.ImposeFields(t,X0,Y0,R,vs);
    Waves.Advection();
    frame = outName+to_string(t)+extension;
    if(t%2==0) Waves.Print(frame.c_str(),Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R, Disk.GetDs(), EpSolver);
    Disk.Fx(Waves); Disk.Fy(Waves);
    F_amp = sqrt(Disk.GetFx()*Disk.GetFx() + Disk.GetFy()*Disk.GetFy());
    if(F_amp > F_max) F_max = F_amp;
    //Disk.UpdatePEFRL(Waves,1.0);
    cout << t << " " << F_max << " " << F_amp << endl;
    //if(t==0)EpSolver.ShowEpsilon();cout << endl;cout << "-------------------" << endl;
    /**if(t > 2*T)
    {
	for(ix=1;ix<Lx-1;ix++)
	{
		pix = Waves.rho(ix,iy,true);	
		if(pix > p_max[ix]) p_max[ix] = pix;
		if(pix < p_min[ix]) p_min[ix] = pix;
	}
    }**/
  }
  
  //Show
	//for(ix=0;ix<Lx;ix++) cout << ix << " " << p_max[ix] << " " << p_min[ix] << endl;
//  Waves.Print("Waves2D.dat",Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R, Disk.GetDs(), EpSolver);
 /**for(int k=0;k<Nodes;k++)
    {
      cout << Disk.GetDotsX()[k] << " " << Disk.GetDotsY()[k] << endl;
    }
  cout << Disk.GetDotsX()[0] << " " << Disk.GetDotsY()[0] << endl;
  **/
  return 0;
}  
