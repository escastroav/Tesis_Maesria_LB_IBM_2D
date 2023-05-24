#include "Waves_D2Q5.h"
#include "IBM.h"
//#include "ComputeEpsilon.h"
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
  /**
  21	22	23	24	25	26

  20	7	8	9	10	27
    
  19	6	1	2	11	28
    
  18	5	0	3	12	29
    
  17	4	15	14	13	30

  16	35	34	33	32	31
  **/
  //first order:
  nb_x[0] = floor(x);	nb_y[0] = floor(y);
  nb_x[1] = nb_x[0];	nb_y[1] = nb_y[0]+1;
  nb_x[2] = nb_x[0]+1;	nb_y[2] = nb_y[0]+1;
  nb_x[3] = nb_x[0]+1;	nb_y[3] = nb_y[0];
  //second order:
  nb_x[4] = floor(x)-1;	nb_y[4] = floor(y)-1;
  //left side:
  nb_x[5] = nb_x[4];	nb_y[5] = nb_y[4]+1;
  nb_x[6] = nb_x[4];	nb_y[6] = nb_y[4]+2;
  nb_x[7] = nb_x[4];	nb_y[7] = nb_y[4]+3;
  //up side:
  nb_x[8] = nb_x[4]+1;	nb_y[8] = nb_y[7];
  nb_x[9] = nb_x[4]+2;	nb_y[9] = nb_y[7];
  nb_x[10] = nb_x[4]+3;	nb_y[10] = nb_y[7];
  //right side:
  nb_x[11] = nb_x[10];	nb_y[11] = nb_y[4]+2;
  nb_x[12] = nb_x[10];	nb_y[12] = nb_y[4]+1;
  nb_x[13] = nb_x[10];	nb_y[13] = nb_y[4];
  //down side: 
  nb_x[14] = nb_x[4]+2;	nb_y[14] = nb_y[4];
  nb_x[15] = nb_x[4]+1;	nb_y[15] = nb_y[4];
  /**third order:
  nb_x[16] = floor(x)-2;nb_y[16] = floor(y)-2;
  //left side:
  nb_x[17] = nb_x[16];	nb_y[17] = nb_y[16]+1;
  nb_x[18] = nb_x[16];	nb_y[18] = nb_y[17]+2;
  nb_x[19] = nb_x[16];	nb_y[19] = nb_y[18]+3;
  nb_x[20] = nb_x[16];	nb_y[20] = nb_y[19]+4;
  nb_x[21] = nb_x[16];	nb_y[21] = nb_y[20]+5;
  //up side:
  nb_x[22] = nb_x[16]+1;nb_y[22] = nb_y[21];
  nb_x[23] = nb_x[16]+2;nb_y[23] = nb_y[21];
  nb_x[24] = nb_x[16]+3;nb_y[24] = nb_y[21];
  nb_x[25] = nb_x[16]+4;nb_y[25] = nb_y[21];
  nb_x[26] = nb_x[16]+5;nb_y[26] = nb_y[21];
  //right side:
  nb_x[27] = nb_x[26];	nb_y[27] = nb_y[16]+4;
  nb_x[28] = nb_x[26];	nb_y[28] = nb_y[16]+3;
  nb_x[29] = nb_x[26];	nb_y[29] = nb_y[16]+2;
  nb_x[30] = nb_x[26];	nb_y[30] = nb_y[16]+1;
  nb_x[31] = nb_x[26];	nb_y[31] = nb_y[16];
  //down side: 
  nb_x[32] = nb_x[16]+4;nb_y[32] = nb_y[16];
  nb_x[33] = nb_x[16]+3;nb_y[33] = nb_y[16];
  nb_x[34] = nb_x[16]+2;nb_y[34] = nb_y[16];
  nb_x[35] = nb_x[16]+1;nb_y[35] = nb_y[16];
  **/
}


//--------------- Global Functions ------------

int main(int argc, char * argv[]){
  LatticeBoltzmann Waves;
  double R = 10.0; // 1 dx => 3 um, 1 dt = 2 ns
  double ds = 1.0;  
  int Nodes = (int)( 2 * M_PI * R / ds );
  double Volume = M_PI * R * R; 
  double vs = 0.559016994374948;//(1700/1500) *C; //polystyrene / water with C_o = 1/3
  double density = 500 / Volume;//(1050/1000) / C2 ; // polystyrene / water with rho = B / C^2
  double B = density * vs * vs; // in the code is defined B_o = 1 for water
  double mass = density * Volume;
  int X0 = Lx/4;
  int Y0 = Ly/2;
  double Po = 1/C2;
  IBMDisk Disk(Nodes, R, B, mass, vs, X0, Y0);
  int t,tmax=10*T;
  double k = 2 * M_PI / Lambda;
  double F_lin=0, F_lin_max=0, F_lin_min=0;
  double F_sqr=0, F_sqr_max=0, F_sqr_min=0;
  double rho0=Po,Jx0=0,Jy0=0;
  string outName = "Wav3D_st=";
  string extension = ".dat";
  string frame = "";
  int iy = Ly/2; int ix = 0;
    
  cout.precision(8);
  //cout << "P_o" << "\t" << "R" << "\t" << "k" << "\t" << "F_min" << "\t" << "F_max" << "\t" << "<F>" << endl;
  //Start
  Waves.Start(rho0,Jx0,Jy0,0,0,X0,Y0,R,vs);
  //Run
  for(t=0;t<tmax;t++){
    Waves.Collision(Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R, Disk.GetDs(), vs, t);
    Waves.ImposeFields(t,X0,Y0,R,vs,Po);
    Waves.Advection();
    //X positions are relative to the wavelength!
    //Disk.Fx(Waves); 
    F_lin = Disk.Fx_p(Waves);
    F_sqr = Disk.Fx_p2_v2(Waves) + Disk.Fx_J(Waves);
    if(F_lin > F_lin_max && t > tmax - T - 1) F_lin_max = F_lin;
    if(F_lin < F_lin_min && t > tmax - T - 1) F_lin_min = F_lin;
    if(F_sqr > F_sqr_max && t > tmax - T - 1) F_sqr_max = F_sqr;
    if(F_sqr < F_sqr_min && t > tmax - T - 1) F_sqr_min = F_sqr;
    //if(FJ_amp > FJ_max && t >= tmax - T - 1) FJ_max = FJ_amp;
    //if(FJ_amp < FJ_min && t >= tmax - T - 1) FJ_min = FJ_amp;
    //if(t >= tmax - T - 1) F_add += F_tot / T;
    //if(t >= tmax - T - 1) FJ_add += F_sqr / T;
    //if(t >= tmax - T - 1) FL_add += F_lin / T;
    //if(t >= 3*T) FJ_add += FJ_amp / T;
    //Disk.MeasureForce(Waves, t, tmax, T, F_min, F_max, F_add);
    /**cout << scientific
         << t/T << "\t"
	 << F_lin << "\t"
	 << F_sqr << "\t"
	 << 0.5*(F_lin_max + F_lin_min) << "\t"
	 << 0.5*(F_sqr_max + F_sqr_min) << endl;**/
    //Disk.UpdatePEFRL(Waves,1.0);
    //frame = outName+to_string(t)+extension;
    //if(t%2==0) Waves.Print(frame.c_str(),Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R, Disk.GetDs(),vs);
  }
  
  //Show
  //cout << scientific
    //<< Po << "\t"
       //<< R << "\t"
       //<< k << "\t"
       //<< F_min << "\t" 
       //<< F_max << "\t"
       //<< F_add << endl;// << "\t"
       //<< FJ_add << endl;
  cout << scientific
       << vs << "\t"
       << 0.5*(F_lin_max + F_lin_min) << "\t"
       << 0.5*(F_sqr_max + F_sqr_min) << endl;
  //Waves.PrintBoundary("data_over_surface.dat", Nodes, Disk.GetDotsX(), Disk.GetDotsY(), Disk.GetX(), Disk.GetY());
  return 0;
}  
