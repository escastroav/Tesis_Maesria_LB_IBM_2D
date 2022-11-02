#include "Waves_D2Q5.h"
#include "IBM.h"

//--------------- Global Functions ------------

int main(void){
  LatticeBoltzmann Waves;
  int Nodes = 12;double R = 10.0;int X0 = Lx/2;int Y0 = Ly/2;
  IBMDisk Disk(Nodes, R, X0, Y0);
  int t,tmax=200;
  double rho0=0,Jx0=0,Jy0=0;
  
  //Start
  Waves.Start(rho0,Jx0,Jy0,0,0);
  //Run
  for(t=0;t<tmax;t++){
    //Disk.UpdatePEFRL(Waves,1.0);
    Waves.Collision(Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R);
    if(t<10.0 / C)Waves.ImposeFields(t);
    Waves.Advection();
    cout << Disk.GetX() << " " << Disk.GetY() << endl;
  }
  //Show
  Waves.Print("Waves2D.dat",Nodes,Disk.GetDotsX(), Disk.GetDotsY(),Disk.GetBulk(), Disk.GetX(), Disk.GetY(), Disk.GetVx(), Disk.GetVy(), R);
 
  return 0;
}  
