//g++ -o run square_lattice.cc -I/Users/mohammad/include -DWITH_LAPACK -lblas -llapack -g -O2
#include "math.h"
#include <iostream>
#include <iosfwd>
#include <fstream>
#include <math.h>
#include <dmtk/vector.h>
#include <dmtk/matrix.h>
#include <LAPACK/zgeev.h>

#define PI (3.141592653589793)
#define IBITS(n,i) (((n) & 1 << i) >> i)
#define IBSET(n,i) ((n) | (1<<i))
#define IBCLR(n,i) ((n)^(IBITS(n,i) << i))

using namespace std;
using namespace dmtk;

unsigned factorial(unsigned);
void Diagonalize(int, Matrix<complex<double> >, Matrix<complex<double> > &, Vector<complex<double> > &);

int main()
{
	double t, B;
	int L;

	cout << "L*L Lattice, L = ";
	cin >> L;
	t = 1.;
	
	ofstream gout("energy.txt",std::ios::out);
	int dim = L*L;
	double phi0 = 6.6260693e-34 / 1.60217653e-19; // phi0=h/e 
	
	for (double phi = 0; phi < phi0; phi += phi0/60.) { // phi=B*a^2
		Matrix<complex<double> > H(dim,dim);
		Matrix<complex<double> > EigVec(dim,dim);
		Vector<complex<double> > eigen(dim);
		
		for (int site = 0; site < dim; site++) {
			int i = site/L; // row
			int j = site % L; // column
			if ( (site-L) > -1)
				H(site,site-L) = -t * exp(std::complex<double> (0,-2*PI*double(j)*phi/phi0));
			if ( (site+L) < dim)
				H(site,site+L) = -t * exp(std::complex<double> (0,2*PI*double(j)*phi/phi0));
			if ( j > 0)
				H(site,site-1) = -t;
			if ( j < L-1)
				H(site,site+1) = -t;
		}
		// solve the generalized eigenvalue problem
		Diagonalize(dim,H,EigVec,eigen);
		
		for (int i = 0; i < dim; i++)
			gout << phi/phi0 << "	" << eigen(i).real() << endl;
	}
	return 0;
}

// For more info about zgeev go to http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen.html#ga80bb116b6dacd91ef3241373469dd7f8
//=============================================================================
//diagonalize a symmetric matrix -- eigenvalues and eigenvectors
void Diagonalize(int dim, Matrix<complex<double> > mat, Matrix<complex<double> > &evect, Vector<complex<double> > &eval)
{
	int N; N=dim;
	Vector<complex<double> > work(4*N); int info;
	Vector<complex<double> > vals(N); 
	Vector <double> RWORK(2*N);
	Matrix<complex<double> > mts(N,N);
	Matrix<complex<double> > VL(N,N);
	Matrix<complex<double> > VR(N,N);
	
	for(int i1=0;i1<N;i1++){for(int j1=i1+1;j1<N;j1++){
		if( (abs(mat(j1,i1).real()- mat(i1,j1).real()) > 1.e-10) || (abs(mat(j1,i1).imag() + mat(i1,j1).imag()) > 1.e-10)) cout << "ERROR " << "i1:"<< i1 << "	j1:"<< j1 << "	mat(i1,j1):"<< mat(i1,j1) << " mat(j1,i1):" << mat(j1,i1) << endl ; }}
	
	for(int i1=0;i1<N;i1++){for(int j1=0;j1<N;j1++){
		mts(j1,i1) = mat(j1,i1); }}
	
	//cout << " Hamiltonian Dimension: " << dim << " --- " << mat.rows() << endl;
	
	zgeev_('N','N',N,mts.array(),N,vals.array(), VL.array() ,N, VR.array() ,N, work.array(),work.size(),RWORK.array(), info);
	
	for(int i1=0;i1<N;i1++){for(int j1=0;j1<N;j1++){
		evect(i1,j1) = mts(i1,j1);}
		eval(i1) = vals(i1);}
	//cout << " Diagonalization info " << info << " ************ " << endl;
} //End Diagonalization

// Factorial
unsigned factorial(unsigned n)
{
    if (n == 0)
        return 1;
    else
        return n * factorial(n - 1);
}
