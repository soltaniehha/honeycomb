//g++ -o run.exe mappingSquare.cc -I/Users/mohammad/include -DWITH_LAPACK -lblas -llapack -g -O2

#include <iostream>
#include <iosfwd>
#include <fstream>
#include <iomanip>
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

void Diagonalize_complex(int, Matrix<complex<double> >, Matrix<complex<double> > &, Vector<complex<double> > &);

int main()
{
	double B, phi;
	int L, orbital, o, i, j, k, l;
    complex <double> psi_psi, psi_psi_1;
	cout << "L*L Lattice, *L MUST be odd*  L = ";
	cin >> L;
	double phi0 = 6.6260693e-34 / 1.60217653e-19; // phi0=h/e
	double a = 1.; //a: lattice constant
	double alpha = 0.; //the constant in the coulomb potential
	double t = 1.;
	int dim = L*L;
    int centerSite = (L+1) * int(L/2); // Starting from the 'site' in the middle of the square
    double err = 1.e-20;
    double err2 = 1.e-5;
	
	ofstream gout("b_mappingSquare.txt",std::ios::out);
	ofstream phiout("energy_mappingSquare.txt",std::ios::out);
		
	// Labeling the sites that we use to build complete orbitals in square lattice.
	Vector<int > os(dim); // occupied site
	int mid_point = int(L/2);
	for (int i = 0; i < L; i++) {
		int valid_points, start_point;
		if ( i < mid_point) {
			valid_points = 2 * i + 1;
			start_point = mid_point - i;
		} else {
			valid_points = L - (2 * (i-mid_point));
			start_point = mid_point - (L-i-1);
		}
		
		for (int j = 0; j < valid_points; j++) {
			os(i*L + start_point) = 1;
			start_point++;
		}
	}
	
    // Changing the magnetic field in order to find the change in energy
    for (phi = 0; phi < phi0; phi += phi0/60.) {
        orbital = 0;
        cout << endl << "===================> phi/phi0: " << phi/phi0 << endl << endl;
        Matrix<complex<double> > psi(1,dim); // Will resize later on
        Vector<complex<double> > b2(1);
        Vector<complex<double> > a(1);
        psi(0,centerSite) = exp(std::complex<double> (0,0)); // exp(std::complex<double> (0,0)) = 1. ------ creating the centerSite!
       // o = 0;
        b2(0) = 0.;
        while (true) 
		{ // Building Psi(orbitals)
			Vector<complex<double> > temp_psi(dim);
            for (k = 0; k < dim ; k++)
                if (abs(psi(orbital,k)) > 1.e-200) {
                    i = k/L; // column
                    j= k % L; // row
                    if ( os(k-L) == 1 ) // to the left
                        temp_psi(k-L) += -t * exp(std::complex<double> (0,-2*PI*double(j)*phi/phi0)) * psi(orbital,k);
                    if ( os(k+L) == 1 ) // to the right
                        temp_psi(k+L) += -t * exp(std::complex<double> (0,2*PI*double(j)*phi/phi0)) * psi(orbital,k);
                    if ( os(k-1) == 1 ) // down
                        temp_psi(k-1) += -t * psi(orbital,k);
                    if ( os(k+1) == 1 ) // up
                        temp_psi(k+1) += -t * psi(orbital,k);
                }
            
            psi_psi = 0; // psi_psi = <psi_i|psi_i>
			psi_psi_1 = 0; // psi_psi = <psi_i|psi_i>
            for (j= 0; j < dim; j++)
                psi_psi += psi(orbital,j) * conj(psi(orbital,j));
            
			// calculating 'a'
            complex <double> psi_H_psi = 0; // psi_H_psi = <psi_i|H_band|psi_i>
            for (k= 0; k < dim; k++)
                if (k != centerSite)
                    if (abs(psi(orbital,k)) > 1.e-10) {
                        int i_site = centerSite/L; // column of the center point
                        int j_site = centerSite % L; // row of the center point
                        i = k/L; // column
                        j= k % L; // row
                        double a_i;
                        a_i = alpha / sqrt(double((i_site-i)*(i_site-i) + (j_site-j)*(j_site-j)));
                        psi_H_psi += a_i * psi(orbital,k) * conj(psi(orbital,k));
                    }

            a(orbital) = psi_H_psi/psi_psi; 

			// Calculating 'b'
           	psi_psi_1 = 0;
            psi_psi = 0;
			complex <double> b2_temp;
			if (orbital > 0) {
                for (l= 0; l < dim; l++) {
                    psi_psi_1 += psi(orbital-1,l) * conj(psi(orbital-1,l));
                    psi_psi += psi(orbital,l) * conj(psi(orbital,l));
                }
                b2_temp = psi_psi/psi_psi_1;
			}
            
            if (( (abs(b2_temp) > err) || orbital == 0) && orbital < int(L/2)) {
                gout  << orbital << "	" << sqrt(b2_temp.real()) << endl;
                
                psi.resize(orbital+2,dim);
                b2.resize(orbital+2);
                a.resize(orbital+2);
                
				b2(orbital) = b2_temp;
                for (k = 0; k < dim ; k++)
                    psi(orbital+1,k) += temp_psi(k);
                
                // Adding the two terms with coef a and b to the new psi/orbital.
                
                // Subtracting "b" term
                if (orbital > 0)
                    for (k = 0; k < dim; k++)
                        psi(orbital+1,k) += -b2(orbital) * psi(orbital-1,k);
                
                // Subtracting "a" term
                for (k = 0; k < dim; k++)
                    psi(orbital+1,k) += -a(orbital) * psi(orbital,k);
                
                    cout << "orbital:" << orbital << "  b^2:" << fixed << setprecision(5) << b2(orbital).real() << "    b:" << fixed << setprecision(5) << sqrt(b2(orbital).real()) << "   a:" << fixed << setprecision(5) << a(orbital).real() << endl;

                // ORTHOGONALIZING |ORBITAL+1> WITH THE REST OF ORBITALS

                for (o = 0; o < orbital+1; o++) {
                    complex <double> o_dot_orbital = 0; // < o | orbital+1 >
                    for (k= 0; k < dim; k++)
                        o_dot_orbital += psi(orbital+1,k) * conj(psi(o,k));

                    if (abs(o_dot_orbital) > err) {
                        complex <double> o_dot_o = 0;
                        for (k= 0; k < dim; k++)
                            o_dot_o += psi(o,k) * conj(psi(o,k));
                        
                        for (k = 0; k < dim; k++)
                                psi(orbital+1,k) -= psi(o,k) * o_dot_orbital / o_dot_o;
                        
	                    complex <double> o_dot_orbital_new = 0; // < o | orbital+1 >
                        for (k= 0; k < dim; k++)
                            o_dot_orbital_new += psi(o,k) * conj(psi(orbital+1,k));
						if ( (abs(o_dot_orbital) > err2) && (abs(o_dot_orbital_new) > abs(o_dot_orbital)) )
	                        cout << "^^^^ the overlap between " << orbital+1 << " and " << o <<  " was:" << abs(o_dot_orbital) << ", and now is:" << abs(o_dot_orbital_new) << endl;
                    }
                }

                
                psi_psi = 0;
                for (k= 0; k < dim; k++)
                    psi_psi += psi(orbital,k) * conj(psi(orbital,k));
                // Confirming that psi(orbital+1) is orthogonal to all of the other orbitals
                for (l= 0; l < orbital+1; l++) {
                    psi_psi_1 = 0;
                    for (k= 0; k < dim; k++)
                        psi_psi_1 += psi(l,k) * conj(psi(l,k));
                    complex <double> overlap = 0;
                    for (j= 0; j < dim; j++)
                      	overlap += psi(l,j) * conj(psi(orbital+1,j));
                    overlap = overlap / psi_psi_1 / psi_psi;
                    if (abs(overlap) > err2)
                        cout << "       < " << l << " | " << orbital+1 << " > = " << abs(overlap) << endl;
                }
				orbital++;
			} else
                break;
        }
		
        //Normalizing all Psi s
        for (o = 0; o < orbital; o++) {
            psi_psi = 0;
            for (j= 0; j < dim; j++)
                psi_psi += psi(o,j) * conj(psi(o,j));
            for (j= 0; j < dim; j++)
                psi(o,j) = psi(o,j) / sqrt(psi_psi);
        }
        
        Matrix<complex<double> > H(orbital,orbital);
		Matrix<complex<double> > EigVec(orbital,orbital);
		Vector<complex<double> > eigen(orbital);
		
        for (o = 0; o < orbital; o++) {
            H(o,o) = a(o);
            if (o < orbital-1)
                H(o,o+1) = sqrt(b2(o+1).real());
            if(o > 0)
                H(o,o-1) = sqrt(b2(o).real());
        }
			
		// solve the eigenvalue problem 
		Diagonalize_complex(orbital,H,EigVec,eigen);
		
        cout << endl;
		for (o= 0; o < orbital; o++)
			phiout << phi/phi0 << "	" << eigen(o).real() << endl;
            
		cout << endl << "***********************************************************************" << endl << "***********************************************************************" << endl << endl;
	}
	return 0;
}

//diagonalize a symmetric complex matrix -- gives eigenvalues and eigenvectors using LAPACK library zgeev.
//For more info: http://www-heller.harvard.edu/people/shaw/programs/lapack.html
void Diagonalize_complex(int dim, Matrix<complex<double> > mat, Matrix<complex<double> > &evect, Vector<complex<double> > &eval)
{
	char  jobvl = 'N'; /* V/N to calculate/not calculate the left eigenvectors
		  of the matrix H.*/
	char  jobvr = 'V'; // As above, but for the right eigenvectors.
	int N; N=dim;
	int lwork = 4 * dim;
	Vector<complex<double> > work(lwork);
	int info;
	Vector<complex<double> > vals(dim); 
	Vector <double> RWORK(2*dim);
	Matrix<complex<double> > mts(dim,dim);
	Matrix<complex<double> > VL(dim,dim);
	Matrix<complex<double> > VR(dim,dim);

	for (int i1 = 0; i1 < dim; i1++) {
		for (int j1 = i1 + 1; j1 < dim; j1++) {
			if( (abs(mat(j1,i1).real()- mat(i1,j1).real()) > 1.e-10) || (abs(mat(j1,i1).imag() + mat(i1,j1).imag()) > 1.e-10)) cout << "ERROR " << "i1:"<< i1 << "	j1:"<< j1 << "	mat(i1,j1):"<< mat(i1,j1) << " mat(j1,i1):" << mat(j1,i1) << endl ; 
		}
	}
	
	int i, j, k;
	for (i = 0; i < dim; i++){
		for (j = 0; j < dim; j++){
			mts(j,i) = mat(j,i); 
		}
	}
		
	cout << endl << " Hamiltonian Dimension: " << dim << " --- " << mat.rows() << endl;
		
	zgeev_(jobvl, jobvr, dim, mts.array(), dim, vals.array(), VL.array() ,dim, VR.array() , dim, work.array(), lwork, RWORK.array(), info);
	
	for (i = 0; i < dim; i++) {
		for(j = 0; j < dim; j++)
			evect(i,j) = VR(i,j);
		eval(i) = vals(i);
	}
	
// Sorting the eigenstates and eigenenergies in an increasing order!
		
	complex<double> temp;
	for (j = 0; j < dim; j++)
		for (i = 0; i < dim-1; i++)
			if (eval(i).real() > eval(i+1).real()) {
				temp = eval(i); eval(i) = eval(i+1); 
				eval(i+1) = temp;
				
				for (k = 0; k < dim; k++) {
					temp = evect(i,k);
					evect(i,k) = evect(i+1,k);
					evect(i+1,k) = temp;
				}
			}
		
	cout << " Done! -- Diagonalization info " << info << " ************ " << endl << endl;
	
} //End Diagonalization