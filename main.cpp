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
  double R = 10.0; // 1 dx => 1 um, 1 dt = 1 ns
  double Gamma = 1.0;
  double a0 = 10.0;
  double b0 = 10.0;
  double Phi = M_PI * 0.0625 * 4.0;
  double D = 20.0;
  double ds = 1.0;
  int Nodes = (int)( 2 * M_PI * R / ds );
  double Volume = M_PI * R * R; 
  double vs = 0.25;//(1699/1500) *C; //polystyrene / water with C_o = 1/3
  double density = 4.0;//500 / Volume;//(1050/1000) / C2 ; // polystyrene / water with rho = B / C^2
  double B = density * vs * vs; // in the code is defined B_o = 1 for water
  double mass = density * Volume;
  int X1 = Lx/2 - D*cos(Phi);
  int X2 = Lx/2 + D*cos(Phi);
  int Y1 = Ly/2 - D*sin(Phi);
  int Y2 = Ly/2 + D*sin(Phi);
  double Po = 10.0/C2;
  IBMDisk Disk1(Nodes, R, Gamma, Phi, a0, b0, B, mass, vs, X1, Y1);
  IBMDisk Disk2(Nodes, R, Gamma, Phi, a0, b0, B, mass, vs, X2, Y2);
  LatticeBoltzmann Waves(&Disk1, &Disk2);
  int t,tmax=100 * T;
  double k = 2 * M_PI / Lambda;
  double T_z=0, T_z_add=0, ART_z = 0;
  double F_x_1=0, F_x_1_add=0, ARF1_x = 0;
  double F_y_1=0, F_y_1_add=0, ARF1_y = 0;
  double F_x_2=0, F_x_2_add=0, ARF2_x = 0;
  double F_y_2=0, F_y_2_add=0, ARF2_y = 0;
  double rho0=Po,Jx0=0,Jy0=0;
  string outName = "speedDisk1Disk2_st=";
  string extension = ".dat";
  string frame = "";
  int iy = Ly/2; int ix = 0;
  double dtMD = 1.0;
  cout.precision(8);
  //Start
  Waves.Start(rho0,Jx0,Jy0,0,0);//,X0,Y0,R,vs);
  //Run
  for(t=0;t<tmax;t++){
    Waves.Collision();//(Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), Phi, R, a0, b0, Disk.GetDs(), vs, t);
    Waves.ImposeFields(t,Po);//,X0,Y0,R,vs,Po);
    Waves.Advection();
    //X positions are relative to the wavelength!
    //F_x = Disk.Fx_p2_v2(Waves) + Disk.Fx_J(Waves);
    if(t > 10*T){
    F_x_1 = Disk1.Fx(Waves);
    F_x_1_add+=F_x_1;
    F_y_1 = Disk1.Fy(Waves);
    F_y_1_add+=F_y_1;

    F_x_2 = Disk2.Fx(Waves);
    F_x_2_add+=F_x_2;
    F_y_2 = Disk2.Fy(Waves);
    F_y_2_add+=F_y_2;
    //if(F_x > F_x_max) F_x_max = F_x;
    //if(F_x < F_x_min) F_x_min = F_x;
    //ARF_x = 0.5*(F_x_max + F_x_min);
    //F_y = Disk.Fy_p2_v2(Waves) + Disk.Fy_J(Waves);
    //F_y = Disk.Fy(Waves);
    //if(F_y > F_y_max) F_y_max = F_y;
    //if(F_y < F_y_min) F_y_min = F_y;
    //ARF_y = 0.5*(F_y_max + F_y_min);
    //T_z = Disk.Tz_J(Waves) + Disk.Tz_p2_v2(Waves);
    //T_z = Disk1.Tz(Waves, D, Phi);
    //T_z_add+=T_z;
    //if(T_z > T_z_max) T_z_max = T_z;
    //if(T_z < T_z_min) T_z_min = T_z;
    //ART_z = 0.5*(T_z_max + T_z_min);
    //cout << t/T << " " << F_x << " " << F_x_add / T << " " << ARF_x << endl;
    //cout << t/T << " " << T_z << " " << T_z_add / T << " " << ART_z << endl;
    if(t % (int)T == 0){
    Disk1.UpdatePEFRL(ARF1_x,ARF1_y,dtMD);
    Disk2.UpdatePEFRL(ARF2_x,ARF2_y,dtMD);
    //ARF_x = 0.5*(F_x_max + F_x_min);
    ARF1_x = F_x_1_add / T; 
    ARF2_x = F_x_2_add / T; 
    ARF1_y = F_y_1_add / T; 
    ARF2_y = F_y_2_add / T; 
    //ART1_z = T_z_add / T; 
     cout << t/T << "\t"
          << ARF1_x << "\t"
          << ARF1_y << "\t"
          << ARF2_x << "\t"
          << ARF2_y << "\t"
          << Disk1.GetX() << "\t"
          << Disk1.GetY() << "\t"
          << Disk2.GetX() << "\t"
          << Disk2.GetY() << endl;
    F_x_1_add = 0;
    F_x_2_add = 0;
    F_y_1_add = 0;
    F_y_2_add = 0;
    //T_z_add = 0;
    }}
     //ARF_x = 0.5*(F_x_max + F_x_min);
     //ARF_y = 0.5*(F_y_max + F_y_min);
    //Disk.UpdatePEFRL(ARF_x,ARF_y,dtMD);}//(Waves,dtMD);
    //frame = outName+to_string(t)+extension;
    //if(t%2==0 && t>10*T) Waves.Print(frame.c_str(),Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R, Phi, a0, b0, Disk.GetDs(),vs);
    //Waves.Print(frame.c_str());//,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R, Phi, a0, b0, Disk.GetDs(),vs);
  }
  
  //Show
  /**
  cout << "D" << "\t"
  	<< "F_x" << "\t\t"
      << "T_z" << endl;
  cout << D << "\t"
  	<< ARF_x << "\t"
      << ART_z << endl;**/
  //Waves.PrintBoundary("data_over_surface.dat", Nodes, Disk.GetDotsX(), Disk.GetDotsY(), Disk.GetX(), Disk.GetY());
  return 0;
}  
