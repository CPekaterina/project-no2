#include <src/UnitTest++.h>
#include <iomanip>
#include <iostream>


using namespace std;

void printmatrix(double ** A, int n, int m);

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
            if(i==j)A[i][j]=double(2)/(h*h)+((double(i)+1.)*h)*((double(i)+1.)*h);
            else if(i==j+1 || i==j-1)A[i][j]=-1./(h*h);
            else{
                A[i][j]=0;
            }
        }
    }

    printmatrix(A,n,n);

    return UnitTest::RunAllTests();
}


// for testing purposes: print a matrix
void printmatrix(double ** A, int n, int m)
{

    for(int i=0;i<n;i++)
    {
        cout << "| ";
        for(int j=0; j<m; j++)
        {
            cout << setw(6) << A[i][j] << " ";
        }
        cout << "|" << endl;
    }
}
