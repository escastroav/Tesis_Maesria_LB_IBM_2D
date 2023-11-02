#include "Waves_D2Q5.h"

LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=w[5]=w[6]=1.0/8;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;  Vx[5]=0; Vx[6]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1; Vy[5]=0; Vy[6]=0;
  Vz[0]=0;  Vz[1]=0;  Vz[2]=0;  Vz[3]=0;  Vz[4]=0;  Vz[5]=1; Vz[6]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Lz*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}
double LatticeBoltzmann::rho(int ix,int iy,int iz,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Speed(int ix,int iy,int iz,int X,int Y,int Z,double R,double v)
{
  int ixp = ix - X; int iyp = iy - Y; int izp = iz - Z;
  double r = ixp*ixp + iyp*iyp + izp*izp;
  double w = 1.0/3.0;
  return C - (C-v)*0.5*(1-tanh(w*(r-R*R)));
}
double LatticeBoltzmann::Jx(int ix,int iy,int iz,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jy(int ix,int iy,int iz,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jz(int ix,int iy,int iz,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=Vz[i]*fnew[n0]; else sum+=Vz[i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Pi(int ix,int iy,int iz,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=Vx[i]*Vx[i]*fnew[n0]; else sum+=Vx[i]*Vx[i]*f[n0];
  }
  return sum;
}
// ---------------------- Interpolation -----------------------------------
double LatticeBoltzmann::Interpolate(char field, double x, double y,double z)
{
  int ix, iy, iz;
  //field = "P" for rho, "X" for Jx, "Y" for Jy, else is zero.}
  //double* nb_x = new double[16]; double* nb_y = new double[16]; double* nb_z = new double[16];
  double* nb_x = new double[Q]; double* nb_y = new double[Q]; double* nb_z = new double[Q];
  //Neighbours(nb_x, nb_y, nb_z, x, y, z);
  ix = (int)floor(x); iy = (int)floor(y); iz = (int)floor(z);
  for(int k = 0; k<Q; k++)
  { nb_x[k] = ix+Vx[k]; nb_y[k] = iy+Vy[k]; nb_z[k] = iz+Vz[k]; }	
  double interp = 0;
  switch(field)
    {
    case 'R' :
      for(int k = 0; k < Q; k++)	
	interp += rho(nb_x[k], nb_y[k], nb_z[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y)*Kernel(nb_z[k] - z);
    case 'P' :
      for(int k = 0; k < Q; k++)	
	interp += Pi(nb_x[k], nb_y[k], nb_z[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y)*Kernel(nb_z[k] - z);
      break;
    case 'X' :
      for(int k = 0; k < Q; k++)	
	interp += Jx(nb_x[k], nb_y[k], nb_z[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y)*Kernel(nb_z[k] - z);
      break;
    case 'Y' :
      for(int k = 0; k < Q; k++)	
	interp += Jy(nb_x[k], nb_y[k], nb_z[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y)*Kernel(nb_z[k] - z);
      break;
    case 'Z' :
      for(int k = 0; k < Q; k++)	
	interp += Jz(nb_x[k], nb_y[k], nb_z[k], false) * Kernel(nb_x[k] - x)*Kernel(nb_y[k] - y)*Kernel(nb_z[k] - z);
      break;
    }
  delete [] nb_x; delete [] nb_y; delete [] nb_z;

  return interp;
}
//----------------------equilibrium functions & BGK collision rules------------------------------------
double  LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,double Jz0,int i,double c){
  if(i>0)
    return 4*w[i]*(c*c*rho0+Vx[i]*Jx0+Vy[i]*Jy0+Vz[i]*Jz0);
  else
    return rho0*(1-4*c*c*(1-W0));
}  
void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0,double Jz0,double X, double Y, double Z, double Radius, double c){
  int ix,iy,iz,i,n0; double c0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
        for(i=0;i<Q;i++){ //on each direction
	  n0=n(ix,iy,iz,i);
	  c0 = Speed(ix, iy, iz, X, Y, Z, Radius, c);
	  rho0 *= sin(K*ix);
	  Jx0 *= -C*sin(K*ix);
	  f[n0]=feq(rho0,Jx0,Jy0,Jz0,i,c0);
        }
}
void LatticeBoltzmann::Collision(double X, double Y, double Z, double Radius, double c){
  int ix,iy,iz,i,n0; double rho0,Jx0,Jy0,Jz0,Fx0,Fy0,c0,Rxy2;
  for(ix=0;ix<Lx;ix++) //for each cell
   for(iy=0;iy<Ly;iy++) //for each cell
    for(iz=0;iz<Lz;iz++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,iz,false);
      Jx0= Jx(ix,iy,iz,false);
      Jy0=Jy(ix,iy,iz,false);
      Jz0=Jz(ix,iy,iz,false);
      c0=Speed(ix,iy,iz,X,Y,Z,Radius,c);
      //c0=Speed_Ellipse(ix, iy, X, Y, phi, A0, B0, c);
      //Rxy2 = (ix - X)*(ix - X) + (iy - Y)*(iy - Y);
        //if((ix >= floor(X-Radius-2) && ix <= ceil(X+Radius+2))
	//   &&
	//   (iy >= floor(Y-Radius-2) && iy <= ceil(Y+Radius+2)))
	//if(Rxy2 == Radius)
      /**if(Rxy2 >= (Radius - 2*Ds)*(Radius - 2*Ds) && Rxy2 <= (Radius + 2*Ds)*(Radius+2*Ds))
	  {
	    //Jx0 = Jy0 * (iy - Y) / (ix - X) ; Jy0 = Jx0 * (ix - X) / (iy -Y);
	    //fnew[n(ix,iy,1)] = f[n(ix,iy,3)];
	    //fnew[n(ix,iy,3)] = f[n(ix,iy,1)];
	    //fnew[n(ix,iy,2)] = f[n(ix,iy,4)];
	    //fnew[n(ix,iy,4)] = f[n(ix,iy,2)];
	    Fx0=Fbpx(Ndots, ix, iy, dotsx, dotsy, bulk, Ux, Ds);	 
	    Fy0=Fbpy(Ndots, ix, iy, dotsx, dotsy, bulk, Uy, Ds);	 
	  }
      else{Fx0 = Fy0 = 0;}
      Jx0+=0.5*Fx0; Jy0+=0.5*Fy0;
      **/
      if(ix==Lx-1 || ix==0)
	  {
	    fnew[n(ix,iy,iz,1)] = f[n(ix,iy,iz,3)];
	    fnew[n(ix,iy,iz,2)] = f[n(ix,iy,iz,4)];
	    fnew[n(ix,iy,iz,3)] = f[n(ix,iy,iz,1)];
	    fnew[n(ix,iy,iz,4)] = f[n(ix,iy,iz,2)];
	  }
      else{
        for(i=0;i<Q;i++){ //for each velocity vector
	        n0=n(ix,iy,iz,i);
	        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,Jz0,i,c0);
        }}
    }  
}
void LatticeBoltzmann::ImposeFields(int t,double X, double Y, double Z, double Radius, double c, double Po){
  int i,ix,iy,iz,n0;
  double rho0,Jx0,Jy0,Jz0,c0; 
  //an oscillating source in the middle
  ix=0; iy=0; iz=0;
  for(iy=0;iy<Ly;iy++)
  for(iz=0;iz<Lz;iz++){
  ix=1;	  
  rho0= Po*sin(Omega*t); 
  c0 = Speed(ix,iy,iz,X,Y,Z,Radius,c);
  Jx0= Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,Jz0,i,c0);
    }
  /**ix=Lx-2;	  
  rho0= -Po*sin(Omega*t); 
  c0 = Speed(ix,iy,X,Y,Radius,c);
  Jx0= Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,i,c0);
    }**/
  }
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,i,ixnext,iynext,iznext,n0,n0next;
  double D = 1.0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
     for(iz=0;iz<Lz;iz++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext= (ix+Vx[i]+Lx)%Lx;
	iynext= (iy+Vy[i]+Ly)%Ly;
	iznext= (iz+Vz[i]+Lz)%Lz;
	n0=n(ix,iy,iz,i); n0next=n(ixnext,iynext,iznext,i);
	f[n0next]=fnew[n0];
	//Bounce back with coefficient D
	/**
	fnew[n(0,iy,1)] = D*f[n(0,iy,3)];
	fnew[n(0,iy,2)] = D*f[n(0,iy,4)];
	fnew[n(0,iy,3)] = D*f[n(0,iy,1)];
	fnew[n(0,iy,4)] = D*f[n(0,iy,2)];
	fnew[n(Lx-1,iy,1)] = D*f[n(Lx-1,iy,3)];
	fnew[n(Lx-1,iy,2)] = D*f[n(Lx-1,iy,4)];
	fnew[n(Lx-1,iy,3)] = D*f[n(Lx-1,iy,1)];
	fnew[n(Lx-1,iy,4)] = D*f[n(Lx-1,iy,2)];**/
      }
}
void LatticeBoltzmann::PrintBoundary(const char * NameFile, int Ndots, double * dotsx, double * dotsy, double X, double Y)
{
  ofstream MyFile(NameFile); double nx, ny, nz, rho0, Jx0, Jy0, Jz0, r0;
  for(int k = 0; k < Ndots; k++)
  {
      	nx = dotsx[k];	ny = dotsy[k]; 
        rho0=Interpolate('P', nx, ny, nz);
      	Jx0=Interpolate('X', nx, ny, nz); 	
	Jy0=Interpolate('Y', nx, ny, nz);	
	Jz0=Interpolate('Z', nx, ny, nz);	
        r0=sqrt((nx-X)*(nx-X)+(ny-Y)*(ny-Y));
        MyFile<<scientific<<k<<" "<<rho0<<" "<<Jx0<<" "<<Jy0<<" "<<(Jx0*(nx-X)/r0 + Jy0*(ny-X)/r0)<<endl;
  }
  MyFile.close();
}
void LatticeBoltzmann::Print(const char * NameFile, double X, double Y, double Z, double Radius,double v){
  ofstream MyFile(NameFile); double rho0, Fx0, Fy0, Jx0, Jy0, Jz0, c2; int ix,iy,iz=Lz/2;
  MyFile.precision(8);
   for(iy=0;iy<Ly;iy++){
    for(ix=1;ix<Lx;ix++){
      rho0=rho(ix,iy,iz,true);
      //Jx0=Jx(ix,iy,true); Jy0=Jy(ix,iy,true); Jy0=Jy(ix,iy,true);
      c2 = Speed(ix, iy, iz, X, Y, Z, Radius, v);
      MyFile<<scientific<<rho0; 
      if(ix<Lx-1)MyFile<<" ";
    }
    MyFile<<endl;
  }
  MyFile.close();
}
