#ifndef PHASESPACE_INTEGRALS_H
#define PHASESPACE_INTEGRALS_H

#include "parameters.h"
#include "Arsenal.h"
#include "Matrix_elements_sq.h"
#include "chemical_potential.h"
#include "Formfactor.h"

int Calculate_emissionrates(double* m, Chemical_potential* chempotential_ptr, int channel, string filename);
double set_gausspoints(double* m, double* s_pt, double* s_weight, double** t_pt, double** t_weight, double* E1_pt_standard, double* E1_weight_standard, double** E2_pt, double** E2_weight);

double Integrate_E1(double* m, double* mu, double Eq, double T, int channel, double s, double t, double* E1_pt_standard, double* E1_weight_standard, double** E2_pt_standard, double** E2_weight_standard, double* results);
double Integrate_E2(double* m, double* mu, double Eq, double T, int channel, double s, double t, double E1, double** E2_pt_standard, double** E2_weight_standard, double* results);

inline double viscous_integrand(double s, double t, double E1, double E2, double Eq, double T, double* m, double f0_E1, double f0_E2, double f0_E3);
inline double Bose_distribution(double E, double T, double mu);
inline double deltaf_chi(double p);

#endif
