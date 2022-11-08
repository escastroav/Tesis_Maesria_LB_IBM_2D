#include "ComputeEpsilon.h"

ComputeEpsilon::ComputeEpsilon(int N)
{
  Nodes = N;
  A = Eigen::MatrixXd::Zero(N,N);
  epsilon_sol.resize(N);
  
}

ComputeEpsilon::~ComputeEpsilon()
{
  
}

void ComputeEpsilon::BuildMatrix(double * dotsx, double * dotsy, double ds)
{
  // double* nbxi = new double[16]; double* nbyi = new double[16];
  // double* nbxj = new double[16]; double* nbyj = new double[16];

  for(int i=0; i<Nodes; i++)
    for(int j=i; j<Nodes; j++)
      {
	// Neighbours(nbxi, nbyi, dotsx[i], dotsy[i]);
	// Neighbours(nbxj, nbyj, dotsx[j], dotsy[j]);	
	// for(int k = 0;k<16;k++)
	//   for(int l = 0;l<16;l++)
	//     A(i,j) += Kernel(nbxi[k]-dotsx[i])*Kernel(nbyi[k]-dotsy[i])*Kernel(nbxj[l]-dotsx[j])*Kernel(nbyj[l]-dotsy[j])*ds;
	A(i,j) = Kernel(dotsx[i]-dotsx[j])*Kernel(dotsy[i]-dotsy[j]);
	if(i!=j)A(j,i)=A(i,j);
      }
  // delete [] nbxi; delete [] nbyi;
  // delete [] nbxj; delete [] nbyj;
}

void ComputeEpsilon::SolveA()
{
  Eigen::VectorXd b = Eigen::VectorXd::Ones(Nodes);
  // ShowA();
  // cout << endl;
  // cout << endl;
  epsilon_sol = A.colPivHouseholderQr().solve(b);
  // ShowEpsilon();
  // cout << endl; cout << endl;
  A = Eigen::MatrixXd::Zero(Nodes,Nodes);
}

