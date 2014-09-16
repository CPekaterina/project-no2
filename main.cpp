#include <src/UnitTest++.h>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;

void printmatrix(double ** A, int n, int m);
void jacobi (double **A, double **R, int n);
void jacobi_rot (double s, double c,int k,int l,int n, double **A);
bool max_nondig(double ** A, int n, int* k, int* l);


TEST(WillFail) {
    CHECK(false);
}

int main()
{

    //read in the number n of steps
    int n_step;
    cout << "How many steps should be taken? n_step= " ;
    cin >> n_step;
    cout << endl;

    //size of the matrix is then:

    int n = n_step - 1;
    
    //read in the maximum radius (dimensionless) rho

    double rho_max;
    cout << "Up to what dimensionless radius should the solution be computed? rho_max= " ;
    cin >> rho_max;
    cout << endl;

    //compute step h

    double h = rho_max/double(n_step);

    // write the matrix for Schrödingers equation:

    double **A;
    A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];
    for (int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i==j)A[i][j]=double(2)/(h*h)+((double(i)+1.)*h)*((double(i)+1.)*h); // d_i=2/h^2+rho_i^2 with rho_i=(i+1)*h
            else if(i==j+1 || i==j-1)A[i][j]=-1./(h*h);                            // -e_i=-1/h^2
            else{
                A[i][j]=0;                                                         // nondiagonal elements
            }
        }
    }

    //setup eigenvector-matrix R

    double **R;
    R = new double* [n];
    for (int i = 0; i < n; i++)
        R[i] = new double[n];
    for (int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i==j)R[i][j]=1;
            else R[i][j]=0;
        }
    }


//    printmatrix(A,n,n);
     jacobi(A,R,n);
//    printmatrix(A,n,n);





    return UnitTest::RunAllTests();
}


// for testing purposes: print a matrix
void printmatrix(double ** A, int n, int m)
{

    for(int i=0;i<n;i++)
    {
        cout << "| ";
        for(int j=0; j<=i; j++)
        {
            cout << setw(6) << A[i][j] << " ";
        }
        cout << "|" << endl;
    }
}

//jacobi algorithm

void jacobi (double **A, double **R, int n)
{
    int k=1,l=0,z=0;                        //set k,l to the first lower-triangle matrix element, z count the number of rotations
    while(max_nondig(A,n,&k,&l)==false)     // while the max^2 of one element of lower-triangle is larger than epsilon
    {
        double t1, t2, t, tau, c, s;
        tau=(A[l][l]-A[k][k])/(2.*A[k][l]); // formulas from the lecture notes to obtain cos() and sin()
        t1= -tau+sqrt(1.+tau*tau);
        t2= -tau-sqrt(1.+tau*tau);
        if(abs(t1)>abs(t2)){t=t2;}          //use the smaller of the roots
        else{t=t1;}
        c=1./sqrt(1+t*t);
        s=c*t;
        z++;                                //increase number of iterations

        jacobi_rot(s,c,k,l,n,A);            //rotate A with c and s to put A[k][l] to zero
    }


    //print the eigenvalues and the number of iterations

    cout << "The eigenvalues are: " << endl;
    for (int i =0; i<n; i++)
    {
        cout << A[i][i] << endl;
    }
    cout << "Number of iterations: " << z << endl;
}

//rotation function

void jacobi_rot (double s, double c,int k,int l,int n, double **A)
{
    for(int i=0;i<n;i++) //see formulas from the lecture notes
    {
      if(i!=k && i!=l)
      {
          double tempik = A[i][k];
          double tempil = A[i][l];
          A[i][k]=tempik*c-tempil*s;
          A[i][l]=tempil*c+tempik*s;
      }
     }
      double tempkk=A[k][k];
      double templl=A[l][l];
      double tempkl=A[k][l];
      A[k][k]=tempkk*c*c-2.*tempkl*c*s+templl*s*s;
      A[l][l]=templl*c*c+2*tempkl*c*s+tempkk*s*s;
      A[k][l]=0.0; //hard coding of the result
}


// finds the maximum of the non-diagonal matrix elements in the lower tri-diagonal and tests if it is larger than the tolerance E

bool max_nondig(double ** A, int n, int* k, int* l)
{
   double E = 1.0e-8;                   //tolerance E
   double max = abs(A[*k][*l]);

   for(int i=1;i<n;i++)                    //iteration over the non-diagonal matrix elements in the lower tri-diagonal
   {
       for(int j=0;j<i;j++)
       {
           if(abs(A[i][j])>max)
           {
               max = abs(A[i][j]);
               *k = i;                      // write row of the new maximum to k;
               *l = j;                      // write column of the new maximum to l;
           }
         if(abs(A[i][j])<=E)A[i][j]=0.0; //set values below E to 0
       }
   }

   return((max*max)<E);
}
