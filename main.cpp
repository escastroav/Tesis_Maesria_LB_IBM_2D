#include "Rotor.h"
#include "Waves_D2Q5.h"
//#include "IBM.h"
//#include "ComputeEpsilon.h"
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
  double R1 = 05.0;
  double R2 = 08.0;
  double Gamma1 = 10.0;
  double Gamma2 = 10.0;
  double a0 = 10.0;
  double b0 = 10.0;
  double Phi = M_PI * 0.0625 * 0.0;
  double D1 = 20.0;
  double D2 = 20.0;
  double D = 0.0;
  double ds = 1.0;
  //int Nodes1 = (int)( 2 * M_PI * R1 / ds );
  //int Nodes2 = (int)( 2 * M_PI * R2 / ds );
  double Volume1 = 4.0 * M_PI * R1 * R1 * R1 / 3.0; 
  double Volume2 = 4.0 * M_PI * R2 * R2 * R2 / 3.0; 
  double vs = 0.25;
  double mu0 = 2.0;
  double density = 1.0;
  //double B = density1 * vs * vs;
  double mass1 = density * Volume1;
  double mass2 = density * Volume2;
  int X = Lx/2 - D*cos(Phi);
  int Y = Ly/2 - D*sin(Phi);
  double Rho0 = 1.0;
  double Po = Rho0/C2;
  Rotor rotor(D1, D2, R1, R2, X, Y, 0.0, 0.0, Phi, 0.0, Rho0, mu0, vs);   
  int t,tmax=1000; //*T;
  double k = 2 * M_PI / Lambda;
  double dtMD = 1e+0;
  double OmegaB = M_PI / 10000.0;
  double eta0 = 0.00; 
  double B0 = 10000;
  cout.precision(8);
  //Start
  //Waves.Start(rho0,Jx0,Jy0,0,0);
  //Run
  for(t=0;t<tmax;t++)
  {
   cout << t/dtMD << " "  
	<< rotor.GetX1() << " "  
	<< rotor.GetY1() << " "  
	<< rotor.GetX2() << " "  
	<< rotor.GetY2() << " "  
	<< rotor.GetTheta() << " "  
	<< rotor.GetX() << " "  
	<< rotor.GetY() << endl;  
   rotor.UpdatePEFRL(t/dtMD, dtMD, OmegaB, C, Rho0, k, Po, eta0, B0);
  }
  return 0;
}   
