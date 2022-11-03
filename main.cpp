#include "Waves_D2Q5.h"
#include "IBM.h"

//--------------- Global Functions ------------

int main(void){
  LatticeBoltzmann Waves;
  int Nodes = 150;
  double R = 15.0;
  double mass = 300;
  double vs = 0.1;//C=0.5
  double density = mass / (M_PI*R*R);
  double B = vs*vs*density;
  int X0 = Lx/3;
  int Y0 = Ly/3;
  IBMDisk Disk(Nodes, R, B, mass, X0, Y0);
  int t,tmax=300;
  double rho0=0,Jx0=0,Jy0=0;
  
  //Rotores magnéticos moviéndose en paredes acústicas: Un estudio computacional.
  
  //Start
  Waves.Start(rho0,Jx0,Jy0,0,0);
  //Run
  for(t=0;t<tmax;t++){
    Disk.UpdatePEFRL(Waves,1.0);
    Waves.Collision(Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R);
    Waves.ImposeFields(t);
    Waves.Advection();
    //cout << Disk.GetX() << " " << Disk.GetY() << endl;
  }
  //Show
  Waves.Print("Waves2D.dat",Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R);
  for(int k=0;k<Nodes;k++)
    {
      cout << Disk.GetDotsX()[k] << " " << Disk.GetDotsY()[k] << endl;
    }
  cout << Disk.GetDotsX()[0] << " " << Disk.GetDotsY()[0] << endl;
  return 0;
}  
