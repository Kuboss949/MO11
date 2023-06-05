//
// Created by Kuboss on 26.05.2023.
//

#ifndef MO11_GAUSS_H
#define MO11_GAUSS_H
#include <iostream>
#include <fstream>
#include<math.h>

using namespace std;

void wypiszWektor(double* v, int size);

double * residuum(double* v, double** A, double* b, int n);

double* errorEstimator(const double* prev, const double* next, int n);

bool checkEstimator(double* e);

bool checkResiduum(double* e);

void gaussSeidel(int n, double** A, double *b, double *xo, int iterations);

double* residuumTridiagonal(const double* v, const double* upperDiagonal, const double* lowerDiagonal, const double* diagonal, const double* b, int n);
void gaussSeidelTridiagonal(int n, const double* upperDiagonal, const double* lowerDiagonal, const double* diagonal, double* b, double* xo, int iterations);

#endif //MO11_GAUSS_H
