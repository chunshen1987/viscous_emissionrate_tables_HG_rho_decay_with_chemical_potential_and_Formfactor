#ifndef Physicalconstants_H
#define Physicalconstants_H

#include <cmath>

const double hbarC = 0.197327053;  //GeV*fm
// particle mass
const double mpion = 0.138;
const double mK = 0.494;
const double mrho = 0.770;
const double momega = 0.782;
const double mKstar = 0.892;
const double ma1 = 1.26;
const double Gammaa1 = 0.4;  //decay width for a1

// parameters for massive Yang-Mill theory
const double gamma_parameter = -0.2913;
const double g_tilde =  6.4483;
const double zeta =  0.0585;
const double m0 = 0.875;
const double F_pi_tilde = 0.135;
const double mV = (m0/sqrt(1-gamma_parameter));
const double Z = (sqrt(1-pow(g_tilde*F_pi_tilde,2)/(4*mV*mV)));
const double delta = (1.0 - Z*Z -2*pow(Z,4)*zeta*g_tilde/((1-Z*Z)*sqrt(1-gamma_parameter)));

const double eta1 = ((g_tilde*F_pi_tilde/(2*mV*mV))*sqrt((1-gamma_parameter)/(1+gamma_parameter)) + 4*zeta*Z*Z/(F_pi_tilde*sqrt(1+gamma_parameter)));
const double eta2 = ((g_tilde*F_pi_tilde/(2*mV*mV))*sqrt((1+gamma_parameter)/(1-gamma_parameter)) - 4*gamma_parameter/(F_pi_tilde*g_tilde*sqrt(1-pow(gamma_parameter,2))));
const double C4 = ((pow(g_tilde*F_pi_tilde,2)/(16*pow(mV,4)))*((1+gamma_parameter)/(1-gamma_parameter)) - gamma_parameter/(mV*mV*(1-gamma_parameter))+ 2*gamma_parameter/(pow(g_tilde*F_pi_tilde,2)*(1-gamma_parameter)));

#endif
