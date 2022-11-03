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
double LatticeBoltzmann::Fbpx(int Ndots, int ix, int iy,  double * dotsx, double * dotsy, double bulk, double Ux)
{
  double I_Jx = 0;
  double BU = bulk * Ux;
  double Fx_k = 0;
  double f = 0;
  for(int k=0; k<Ndots; k++){
    I_Jx = Interpolate('X', dotsx[k], dotsy[k]);
    Fx_k = (BU - I_Jx);
    f += Fx_k * Kernel(ix - dotsx[k]) * Kernel(iy - dotsy[k])*1.19;
  }
  return f;
}
double LatticeBoltzmann::Fbpy(int Ndots, int ix, int iy, double * dotsx, double * dotsy, double bulk, double Uy)
{
  double I_Jy = 0;
  double BU = bulk * Uy;
  double Fy_k = 0;
  double f = 0;
  for(int k=0; k<Ndots; k++){
    I_Jy = Interpolate('Y', dotsx[k], dotsy[k]);
    Fy_k = (BU - I_Jy);
    f += Fy_k * Kernel(ix - dotsx[k]) * Kernel(iy - dotsy[k])*1.19;
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
double LatticeBoltzmann::Kernel(double r)
{
  // r = x - x1 / (x0 - x1)
  if(abs(r) < 0.5)
    return (1.0+sqrt(1.0-3.0*r*r))/3.0;
  else if(abs(r) >= 0.5 && abs(r) <= 1.5)
    return (5.0-3.0*abs(r)-sqrt(1.0 - 3.0*(1.0-abs(r))*(1.0-abs(r)) ))/6.0;
  else
    return 0;
}
void LatticeBoltzmann::Neighbours(double* nb_x, double* nb_y, double x, double y)
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
void LatticeBoltzmann::Collision(int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux, double Uy, double Radius){
  int ix,iy,i,n0; double rho0,Jx0,Jy0,Fx0,Fy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false);
      if((ix >= floor(X-Radius-2) && ix <= ceil(X+Radius+2))
	 &&
	 (iy >= floor(Y-Radius-2) && iy <= ceil(Y+Radius+2)))
	{
	  Fx0=Fbpx(Ndots, ix, iy, dotsx, dotsy, bulk, Ux);//cout << Fx0 << " ";
	  Fy0=Fbpy(Ndots, ix, iy, dotsx, dotsy, bulk, Ux);//cout << Fy0 << " ";
	  //cout << endl;
	}
      else{Fx0 = Fy0 = 0;}
      Jx0=Jx(ix,iy,false)+0.5*Fx0; Jy0=Jy(ix,iy,false)+0.5*Fy0;
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
  rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,i);
  }
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
void LatticeBoltzmann::Print(const char * NameFile,int Ndots, double * dotsx, double * dotsy, double bulk, double X, double Y, double Ux, double Uy, double Radius){
  ofstream MyFile(NameFile); double rho0, Fx0, Fy0, Jx0, Jy0; int ix,iy;
  for(ix=1;ix<Lx-1;ix++){
    for(iy=1;iy<Ly-1;iy++){
      rho0=rho(ix,iy,true);
      if((ix >= floor(X-Radius-2) && ix <= ceil(X+Radius+2))
	 &&
	 (iy >= floor(Y-Radius-2) && iy <= ceil(Y+Radius+2)))
	{
	  Fbpx(Ndots, ix, iy, dotsx, dotsy, bulk, Ux);
	  Fbpy(Ndots, ix, iy, dotsx, dotsy, bulk, Ux);
	}
      else{Fx0 = Fy0 = 0;}
      Jx0=Jx(ix,iy,true)+0.5*Fx0; Jy0=Jy(ix,iy,true)+0.5*Fy0;
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
