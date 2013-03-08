#ifndef MATRIX_ELEMENTS_SQ_H
#define MATRIX_ELEMENTS_SQ_H

#include "parameters.h"
#include "Arsenal.h"

double C1p1(double s, double t);
double C1p2(double s, double t);
double C1p3(double s, double t);
double C1p4(double s, double t);
double C1p5(double s, double t);
double C1p6(double s, double t);
double C2p1(double s, double t);
double C2p2(double s, double t);
double C3p1(double s, double t);
double C3p2(double s, double t);
double C4p1(double s, double t);
double C4p2(double s, double t);
double C5p1(double s, double t);
double C5p2(double s, double t);
double C6p1(double s, double t);
double C6p2(double s, double t);
double C7p1(double s, double t);
double C7p2(double s, double t);

double Matrix_elements_sq(double s, double t, int channel);
double Matrix_elements_sq_C1(double s, double t);
double Matrix_elements_sq_C1_omega(double s, double t);
double Matrix_elements_sq_C2(double s, double t);
double Matrix_elements_sq_C3(double s, double t);
double Matrix_elements_sq_C4(double s, double t);
double Matrix_elements_sq_C5(double s, double t);
double Matrix_elements_sq_C6(double s, double t);
double Matrix_elements_sq_C7(double s, double t);

inline double Power(double x, int a);
#endif
