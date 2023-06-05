//
// Created by Kuboss on 26.05.2023.
//

#include "gauss.h"
void wypiszWektor(double* v, int size){
    for (int i = 0; i < size; ++i) {
        cout << v[i] << " ";
    }
}

double * residuum(double* v, double** A, double* b, int n){
    auto *r = new double[n];
    for (int i = 0; i < n; ++i) {
        r[i] = b[i];
        for (int j = 0; j < n; ++j) {
            r[i] -= A[i][j] * v[j];
        }
    }
    return r;
}



bool checkEstimator(double* e){
    if(e[0]<0.00001 && e[1]<0.00001 && e[2]<0.00001 && e[3]<0.00001)
        return true;
    return false;
}

bool checkResiduum(double* e){
    if(e[0]<0.00001 && e[1]<0.00001 && e[2]<0.00001 && e[3]<0.00001)
        return true;
    return false;
}

void gaussSeidel(int n, double** A, double *b, double *xo, int iterations){
    auto* x1=new double[n];
    auto* x=new double[n];
    for(int i=0; i<n; ++i){
        x1[i]=xo[i];
        x[i]=xo[i];
    }
    double* estimator, *res;
    double sum;
    for (int k = 0; k < iterations; ++k) {
        for (int i = 0; i < n; ++i) {
            sum=0;
            for (int j = 0; j < n; ++j) {
                if(i!=j)
                    sum+=A[i][j]*x1[j];
            }
            x1[i]=(b[i]-sum)/A[i][i];
        }

        estimator = errorEstimator(x, x1, n);
        res = residuum(x1, A, b, n);

        /*cout << "Iteracja nr: " << k+1 << "\nPrzyblizenie wyniku: ";
        wypiszWektor(x1, n);
        cout << "\nEstymator bledu rozwiazania: ";
        wypiszWektor(estimator, n);
        cout << "\nResiduum: ";
        wypiszWektor(res, n);
        cout << "\n\n";*/

        for (int i = 0; i < n; ++i) {
            x[i]=x1[i];
        }
        if(k+1==iterations){
            //cout << "Przerwano ze wzgledu na wyczerpanie limitu iteracji.";
        }
        if(checkEstimator(estimator) && checkResiduum(res)){
            //cout << "Przerwano ze wzgledu na kryterium dokladnosci wyznaczenia xn\noraz kryterium wiarygodnosci xn jako przyblizenia pierwiastka\nRozwiazanie: ";
            wypiszWektor(x, n);
            break;
        }
    }


}

double* residuumTridiagonal(const double* v, const double* upperDiagonal, const double* lowerDiagonal, const double* diagonal, const double* b, int n) {
    auto *r = new double[n];
    r[0] = b[0] - diagonal[0] * v[0] - upperDiagonal[0] * v[1];
    for (int i = 1; i < n - 1; ++i) {
        r[i] = b[i] - lowerDiagonal[i-1] * v[i - 1] - diagonal[i] * v[i] - upperDiagonal[i] * v[i + 1];
    }
    r[n - 1] = b[n - 1] - lowerDiagonal[n - 2] * v[n - 2] - diagonal[n - 1] * v[n - 1];
    return r;
}

double* errorEstimator(const double* prev, const double* next, int n){
    auto *e = new double[n];
    for (int i = 0; i < n; ++i) {
        e[i] = fabs(next[i] - prev[i]);
    }
    return e;
}

//funkcja sprawdzająca estymator błędu
bool checkEstimatorNew(const double* prev, const double* next, int n){
    double tol = 0.000001;
    for (int i = 0; i < n; ++i) {
        if(fabs(next[i] - prev[i])>tol){
            return false;
        }
    }
    return true;
}
//funkcja sprawdzająca reziduum
bool checkResiduumNew(const double* v, const double* upperDiagonal, const double* lowerDiagonal, const double* diagonal, const double* b, int n){
    double tol = 0.00000001;
    if((b[0] - diagonal[0] * v[0] - upperDiagonal[0] * v[1])>tol || (b[n - 1] - lowerDiagonal[n - 2] * v[n - 2] - diagonal[n - 1] * v[n - 1])>tol)
        return false;
    for (int i = 1; i < n-1; ++i) {
        if((b[i] - lowerDiagonal[i-1] * v[i - 1] - diagonal[i] * v[i] - upperDiagonal[i] * v[i + 1])>tol){
            return false;
        }
    }
    return true;
}


//funkcja wykonująca algorytm gaussa-seidla na macierzy trójdiagonalnej podanej przez 3 wektory przekątnych
void gaussSeidelTridiagonal(int n, const double* upperDiagonal, const double* lowerDiagonal, const double* diagonal, double* b, double* xo, int iterations) {
    auto* x1 = new double[n];


    for (int k = 0; k < iterations; ++k) {
        for (int i = 0; i < n; ++i) {
            if (i == 0)
                xo[i] = (b[i] - upperDiagonal[i]*x1[i+1])/diagonal[i];
            if (i == (n - 1))
                xo[i] = (b[i] - lowerDiagonal[i-1]*x1[i-1])/diagonal[i];
            else
                xo[i] = (b[i] - (upperDiagonal[i]*x1[i+1]+lowerDiagonal[i-1]*x1[i-1]))/diagonal[i];
        }

        if (k + 1 == iterations) {
            cout << "Nie zbiezne" << '\n';
        }

        if (checkEstimatorNew(xo, x1, n) && checkResiduumNew(x1, upperDiagonal, lowerDiagonal, diagonal, b, n)) {

            for (int i = 0; i < n; ++i) {
                x1[i] = xo[i];
            }
            break;
        }
        for (int i = 0; i < n; ++i) {
            x1[i] = xo[i];
        }

    }
    for (int i = 0; i < n; ++i) {
        b[i] = xo[i];
    }

    delete[] x1;
}
