#include "Waves_D2Q5.h"

LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}
/*void LatticeBoltzmann::SetDisk(IBMDisk & d)
{
  *disk = d;
  }*/
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Fbpx(int Ndots, int ix, int iy,  double * dotsx, double * dotsy, double bulk, double Ux, double ds, ComputeEpsilon & CE)
{
  double I_Jx = 0;
  double BU = bulk * Ux;
  double Fx_k = 0;
  double f = 0;
  //CE.BuildMatrix(dotsx, dotsy,ds); CE.SolveA();
  for(int k=0; k<Ndots; k++){
    I_Jx = Interpolate('X', dotsx[k], dotsy[k]);
    Fx_k = (BU - I_Jx);    
    f += Fx_k * Kernel(ix - dotsx[k]) * Kernel(iy - dotsy[k])*ds*CE.GetEpsilonK(k);
  }
  return f;
}
double LatticeBoltzmann::Fbpy(int Ndots, int ix, int iy, double * dotsx, double * dotsy, double bulk, double Uy, double ds, ComputeEpsilon & CE)
{
  double I_Jy = 0;
  double BU = bulk * Uy;
  double Fy_k = 0;
  double f = 0;
  //CE.BuildMatrix(dotsx, dotsy,ds); CE.SolveA();
  for(int k=0; k<Ndots; k++){
    I_Jy = Interpolate('Y', dotsx[k], dotsy[k]);
    Fy_k = (BU - I_Jy);
    f += Fy_k * Kernel(ix - dotsx[k]) * Kernel(iy - dotsy[k])*ds*CE.GetEpsilonK(k);
  }
  return f;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}
// ---------------------- Interpolation -----------------------------------
double LatticeBoltzmann::Interpolate(char field, double x, double y)
{
  //field = "P" for rho, "X" for Jx, "Y" for Jy, else is zero.}
  double* nb_x = new double[16]; double* nb_y = new double[16];
  Neighbours(nb_x, nb_y, x, y);
  double interp = 0;
  switch(field)
    {
    case 'P' :
      for(int k = 0; k < 16; k++)	
	interp += rho(nb_x[k], nb_y[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y);
      break;
    case 'X' :
      for(int k = 0; k < 16; k++)	
	interp += Jx(nb_x[k], nb_y[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y);
      break;
    case 'Y' :
      for(int k = 0; k < 16; k++)	
	interp += Jy(nb_x[k], nb_y[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y);
      break;
    }
  delete [] nb_x; delete [] nb_y;
  return interp;
}
//----------------------equilibrium functions & BGK collision rules------------------------------------
double  LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i){
  if(i>0)
    return 3*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return rho0*AUX0;
}  
void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0, double Fx0, double Fy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	n0=n(ix,iy,i);
	f[n0]=feq(rho0,Jx0+0.5*Fx0,Jy0+0.5*Fy0,i);
      }
}  
void LatticeBoltzmann::Collision(int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux, double Uy, double Radius, double Ds, ComputeEpsilon & CE){
  int ix,iy,i,n0; double rho0,Jx0,Jy0,Fx0,Fy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false);
      Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      if((ix >= floor(X-Radius-2) && ix <= ceil(X+Radius+2))
	 &&
	 (iy >= floor(Y-Radius-2) && iy <= ceil(Y+Radius+2)))
	{
	  Fx0=Fbpx(Ndots, ix, iy, dotsx, dotsy, bulk, Ux, Ds, CE);//cout << Fx0 << " ";
	  Fy0=Fbpy(Ndots, ix, iy, dotsx, dotsy, bulk, Uy, Ds, CE);//cout << Fy0 << " ";
	}
      else{Fx0 = Fy0 = 0;}
      Jx0+=0.5*Fx0; Jy0+=0.5*Fy0;
      for(i=0;i<Q;i++){ //for each velocity vector
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i);
      }
    }  
}
void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0;
  double lambda,omega,rho0,Jx0,Jy0; lambda=10.0; omega=2*M_PI/lambda*C;
  //an oscillating source in the middle
  ix=1; iy=1;
  for(iy=0;iy<Ly;iy++){
  rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,i);
  }}
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=1;ix<Lx-1;ix++) //for each cell
    for(iy=1;iy<Ly-1;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext= ix+Vx[i]; 
	iynext= iy+Vy[i];
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}
void LatticeBoltzmann::Print(const char * NameFile,int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux, double Uy, double Radius, double ds, ComputeEpsilon & CE){
  ofstream MyFile(NameFile); double rho0, Fx0, Fy0, Jx0, Jy0; int ix,iy;
  for(ix=1;ix<Lx-1;ix++){
    for(iy=1;iy<Ly-1;iy++){
      rho0=rho(ix,iy,true);
      Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      if((ix >= floor(X-Radius-2) && ix <= ceil(X+Radius+2))
	 &&
	 (iy >= floor(Y-Radius-2) && iy <= ceil(Y+Radius+2)))
	{
	  Fx0=Fbpx(Ndots, ix, iy, dotsx, dotsy, bulk, Ux, ds, CE);//cout << Fx0 << " ";
	  Fy0=Fbpy(Ndots, ix, iy, dotsx, dotsy, bulk, Uy, ds, CE);//cout << Fy0 << " ";
	}
      else{Fx0 = Fy0 = 0;}
      Jx0+=0.5*Fx0; Jy0+=0.5*Fy0;
      MyFile<<ix<<" "
	    <<iy<<" "
	    <<rho0<<" "
	    <<Jx0<<" "
	    <<Jy0<<" "
	    <<Fx0<<" "
	    <<Fy0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}
