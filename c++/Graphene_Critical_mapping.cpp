//
//  Graphene_Critical_mapping.cpp
//  
//
//  Created by Andrew Allerdt on 9/24/14.
//
//

#include "Graphene_Critical_mapping.h"

#include <stdio.h>
#include <tgmath.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <valarray>
#include <complex>
#include <iomanip>
#define _USE_MATH_DEFINES
#define PI (3.141592653589793)


using namespace std;

const int l=400;
const int n_states = l;
const int n_sites = (2*l+1)*(2*l+1);
const double eps = 0.0000000000000001;
complex<double> newstates[n_sites][n_states];
complex<double> a[n_states-1];
complex<double> b[n_states-1];

complex<double> VecPrd(complex<double> A[][n_states], int s1, int s2)
{
    complex<double> z(0,0);
    for (int k=0;k<n_sites;k++)
    {
        z = z + conj(A[k][s1])*A[k][s2];
    }
    
    return z;
}

void Normalize(complex<double> B[][n_states], int s1)
{
    complex<double> total(0,0);
    for (int k=0; k<n_sites; k++)
    {
        total = total + B[k][s1]*conj(B[k][s1]);
    }
    total=sqrt(total);
    for (int k=0; k<n_sites; k++)
    {
        B[k][s1]=B[k][s1]/total;
    }
}

void Find(complex<double> B[][n_states],int s1, vector<int>& v2)
{
    for (int i=0;i<n_sites;i++)
    {
        if (abs(B[i][s1])>eps)
        {
            v2.push_back(i);
        }
    }
}

int main()
{
    
//    int l1;
//    double B1;
//    double V1;
//    char *inname = "input";
//    ifstream infile(inname);
//    infile >> l1;
//    infile >> B1;
//    infile >> V1;
//    infile.close();
    
//    cout << "N sites = "<< n_sites << endl;
    int P[n_sites+1][2];
    int seed=n_sites/2 + 1;
//    cout << "Seed = " << seed << endl;
    
    
    double t=1;
    complex<double> t_up;
    complex<double> t_down;
    int a0=1;
    int a_z=1;
    double on_e;
    double dis;
    double B=0.4;
    double V=0.0;
//    cout << "B = " << B << endl;
//    cout << "V = " << V << endl;
    
    double xpos;
    double ypos;
    int active;
    int rm=2*l+1;       // right move
    int lm=-(2*l+1);    // left move
    int um=1;           // up move
    int dm=-1;          // down move
    
    vector<int> occ_sites;
    int hop2sites[3];
    
    
    // Generate Lattice
    int k=1;
    while (k<=n_sites)
    {
        for (int m=-l; m<=l; m++)
        {
            for (int n=-l; n<=l; n++)
            {
                P[k][0]=m;
                P[k][1]=n;
                k++;
            }
        }
    }
    
    newstates[seed][0]=1;

    ypos = sqrt(3)*0.5*a0*P[seed][1];
    xpos = P[seed][0];
    if (active%2==1) // odd site
        xpos=xpos*3*a0/2;
    else if (active%2==1)         // even site
        xpos=(xpos*3*0.5 - 0.5)*a0;
    
    hop2sites[0] = seed+rm; // right
    hop2sites[1] = seed+um; // up neighbor
    hop2sites[2] = seed+dm; // down neighbor
    
    t_up = complex<double>(0,B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
    t_down = complex<double>(0,-B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
    t_up=t*exp(t_up);
    t_down=t*exp(t_down);
    on_e = -V/a_z;
    
    newstates[hop2sites[0]][1] = t*newstates[seed][0] + newstates[hop2sites[0]][1];
    newstates[hop2sites[1]][1] = t_up*newstates[seed][0] + newstates[hop2sites[1]][1];
    newstates[hop2sites[2]][1] = t_down*newstates[seed][0] + newstates[hop2sites[2]][1];
    newstates[seed][1] = on_e*newstates[seed][0] + newstates[seed][1];
    
    a[0] = VecPrd(newstates,0,1);
    
    for (int i=0;i<n_sites;i++)
    {
        newstates[i][1] = newstates[i][1] - a[0]*newstates[i][0];
    }
    Normalize(newstates,1);
    
    // Main Loop
    for (int k=2;k<n_states;k++)
    {
        Find(newstates,k-1,occ_sites);
        
        for (int i=0; i<occ_sites.size(); i++)
        {
            active = occ_sites[i];
            ypos = sqrt(3)*0.5*a0*P[active][1];
            xpos = P[active][0];
            if (active%2==1)
                xpos=xpos*3*a0/2;
            else if (active%2==0)
                xpos=(xpos*3*0.5 - 0.5)*a0;
            
            dis = (xpos*xpos + ypos*ypos);
            on_e = -V/sqrt(a_z*a_z + dis );
            
            if (active%2==1) // odd site
            {
                t_up = complex<double>(0,B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
                t_down = complex<double>(0,-B*(sqrt(3)*a0*xpos/2 - sqrt(3)*a0*a0/8));
                t_up=t*exp(t_up);
                t_down=t*exp(t_down);
                
                hop2sites[0] = active+rm; // right
                hop2sites[1] = active+um; // up neighbor
                hop2sites[2] = active+dm; // down neighbor
                
                newstates[hop2sites[0]][k] = t*newstates[active][k-1] + newstates[hop2sites[0]][k];
                newstates[hop2sites[1]][k] = t_up*newstates[active][k-1] + newstates[hop2sites[1]][k];
                newstates[hop2sites[2]][k] = t_down*newstates[active][k-1] + newstates[hop2sites[2]][k];
                newstates[active][k] = on_e*newstates[active][k-1] + newstates[active][k];
            }
            else if (active%2==0) // even site
            {
                t_up = complex<double>(0,B*(sqrt(3)*a0*xpos/2 + sqrt(3)*a0*a0/8));
                t_down = complex<double>(0,-B*(sqrt(3)*a0*xpos/2 + sqrt(3)*a0*a0/8));
                t_up=t*exp(t_up);
                t_down=t*exp(t_down);
                
                hop2sites[0] = active+lm; // left
                hop2sites[1] = active+um; // up neighbor
                hop2sites[2] = active+dm; // down neighbor
                
                newstates[hop2sites[0]][k] = t*newstates[active][k-1] + newstates[hop2sites[0]][k];
                newstates[hop2sites[1]][k] = t_up*newstates[active][k-1] + newstates[hop2sites[1]][k];
                newstates[hop2sites[2]][k] = t_down*newstates[active][k-1] + newstates[hop2sites[2]][k];
                newstates[active][k] = on_e*newstates[active][k-1] + newstates[active][k];
            }
        }
        
        a[k-1] = VecPrd(newstates,k-1,k);
        b[k-1] = VecPrd(newstates,k-2,k);
        
        for (int m=0;m<n_sites;m++)
        {
            newstates[m][k] = newstates[m][k] - a[k-1]*newstates[m][k-1] - b[k-1]*newstates[m][k-2];
        }
        
        complex<double> overlap;
        for (int j=0;j<k;j++)
        {
            overlap = VecPrd(newstates,j,k);
            if (abs(overlap)>10*eps)
            {
                for (int p=0;p<n_sites;p++)
                    newstates[p][k] = newstates[p][k] - overlap*newstates[p][j];
//                cout << "Overlap=" << overlap << endl;
            }
        }
        
        Normalize(newstates,k);
        occ_sites.clear();
    }
    
    
    
    
    ofstream outfile;
    outfile.open("gch_B=0.4_V=0.txt");
    for (int i=0;i<n_states; i++)
    {
        outfile << setprecision(14) << real(a[i]) << "\t" << setprecision(14) << real(b[i]) << endl;
    }
    outfile.close();
    
//    cout << "B values:"<< endl;
//    for (int i=0;i<n_states;i++)
//        cout << b[i] << endl;
//    cout << " a values:" << endl;
//    for (int i=0;i<n_states;i++)
//        cout << a[i] << endl;
    
    
    return 0;
}


