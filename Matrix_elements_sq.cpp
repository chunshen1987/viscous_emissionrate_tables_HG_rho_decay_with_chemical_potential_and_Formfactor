#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Matrix_elements_sq.h"
using namespace std;

double Matrix_elements_sq(double s, double t, int channel)
{
   double result;
   switch(channel)
   {
      case 1:  //C1: pi + rho -> pi + gamma
         result = Matrix_elements_sq_C1(s, t);
         break;
      case 2:  //C1_omega: pi + rho -> omega -> pi + gamma
         result = Matrix_elements_sq_C1_omega(s, t);
         break;
      case 3:  //C2: pi + pi -> rho + gamma
         result = Matrix_elements_sq_C2(s, t);
         break;
      case 4:  //C4: pi + Kstar -> K + gamma
         result = Matrix_elements_sq_C4(s, t);
         break;
      case 5:  //C5: pi + K -> Kstar + gamma
         result = Matrix_elements_sq_C5(s, t);
         break;
      case 6:  //C6: rho + K -> K + gamma
         result = Matrix_elements_sq_C6(s, t);
         break;
      case 7:  //C7: K + Kstar -> pi + gamma
         result = Matrix_elements_sq_C7(s, t);
         break;
      case 8:  //C3: rho -> pi + pi + gamma
         result = Matrix_elements_sq_C3(s, t);
         break;
      default:
         cout << "Matrix elements squarted ERROR: can not find the corresponding channel! (channel =  " << channel << ")" << endl;
         exit(1);
         
   }
   return(result);
}

double Matrix_elements_sq_C1(double s, double t)
{
  //pi + rho -> pi + gamma
  double isospin_factor_C1p1 = 2.0;
  double isospin_factor_C1p2 = 2.0;
  double isospin_factor_C1p3 = 2.0;
  double result =  (C1p1(s,t)*isospin_factor_C1p1 
                   + C1p2(s,t)*isospin_factor_C1p2 
                   + C1p3(s,t)*isospin_factor_C1p3);
  return (result);  
}

double Matrix_elements_sq_C1_omega(double s, double t)
{
  //pi + rho -> omega -> pi + gamma
  double isospin_factor_C1p4 = 1.0;
  double isospin_factor_C1p5 = 2.0;
  double isospin_factor_C1p6 = 2.0;
  double result = (C1p4(s,t)*isospin_factor_C1p4
                   + C1p5(s,t)*isospin_factor_C1p5
                   + C1p6(s,t)*isospin_factor_C1p6);
  return (result);
}

double Matrix_elements_sq_C2(double s, double t)
{
  //pi + pi -> rho + gamma
  double isospin_factor_C2p1 = 1.0;
  double isospin_factor_C2p2 = 2.0;
  double result = (C2p1(s,t)*isospin_factor_C2p1
                   + C2p2(s,t)*isospin_factor_C2p2);
  return (result);  
}

double Matrix_elements_sq_C3(double s, double t)
{
  //rho -> pi + pi + gamma
  double isospin_factor_C3p1 = 1.0;
  double isospin_factor_C3p2 = 2.0;
  double result = (C3p1(s,t)*isospin_factor_C3p1
                   + C3p2(s,t)*isospin_factor_C3p2);
  return (result);  
}

double Matrix_elements_sq_C4(double s, double t)
{
  //pi + Kstar -> K + gamma
  double isospin_factor_C4p1 = 4.0;
  double isospin_factor_C4p2 = 4.0;
  double result = (C4p1(s,t)*isospin_factor_C4p1
                   + C4p2(s,t)*isospin_factor_C4p2);
  return (result);
}

double Matrix_elements_sq_C5(double s, double t)
{
  //pi + K -> Kstar + gamma
  double isospin_factor_C5p1 = 4.0;
  double isospin_factor_C5p2 = 4.0;
  double result = (C5p1(s,t)*isospin_factor_C5p1
                   + C5p2(s,t)*isospin_factor_C5p2);
  return (result);
}

double Matrix_elements_sq_C6(double s, double t)
{
  //rho + K -> K + gamma
  double isospin_factor_C6p1 = 4.0;
  double isospin_factor_C6p2 = 4.0;
  double result = (C6p1(s,t)*isospin_factor_C6p1
                   + C6p2(s,t)*isospin_factor_C6p2);
  return (result);
}

double Matrix_elements_sq_C7(double s, double t)
{
  //K + Kstar -> pi + gamma
  double isospin_factor_C7p1 = 4.0;
  double isospin_factor_C7p2 = 4.0;
  double result = (C7p1(s,t)*isospin_factor_C7p1
                   + C7p2(s,t)*isospin_factor_C7p2);
  return (result);
}

/*************************************************************************************/
//Matrix elements square
/*************************************************************************************/

/*************************************************************************************/
//pi + rho -> pi + gamma  (C1.1 + C1.2 + C1.3)
/*************************************************************************************/
double C1p1(double s, double t)
{
   double C = 0.059;
   double ghat = g_tilde;
	
   double u = -s - t + 2*mpion*mpion + mrho*mrho; 

   double MaMa = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) + Power(mpion,4)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 4*Power(s,2) - 2*s*t) + 
          Power(s,2)*(-Power(mrho,4) + Power(s,2) - 2*Power(mrho,2)*t + 2*s*t + 2*Power(t,2)) + 
          2*Power(mpion,2)*s*(Power(mrho,4) + Power(mrho,2)*(s + t) - 2*s*(s + t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(Power(mrho,4) + 4*Power(mrho,2)*s + 4*Power(s,2) - 2*s*t) + 
          Power(s,2)*(Power(mrho,4) + Power(s,2) + 2*Power(mrho,2)*(s - t) + 2*s*t + 2*Power(t,2)) - 
          2*Power(mpion,2)*s*(Power(mrho,4) + Power(mrho,2)*(2*s - t) + 2*s*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) - 
          2*Power(mpion,2)*(-Power(mrho,2) + s)*(2*s*(s + t) - Power(mrho,2)*(2*s + t)) + 
          Power(mpion,4)*(3*Power(mrho,4) + 2*s*(2*s - t) + Power(mrho,2)*(-6*s + 2*t)) + 
          s*(-Power(mrho,2) + s)*(Power(s,2) + 2*s*t + 2*Power(t,2) - Power(mrho,2)*(s + 2*t)))))/
   (32.*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

    double MaMb = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (Power(eta2,2)*(Power(mpion,8) - 4*Power(mpion,2)*s*t*(Power(mrho,2) + s + t) + 
          Power(mpion,4)*(2*s*t + Power(mrho,2)*(s + t)) + 
          s*t*(Power(s,2) + 3*s*t + Power(t,2) + Power(mrho,2)*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) + Power(mpion,4)*(2*Power(mrho,4) + 2*s*t - 3*Power(mrho,2)*(s + t)) + 
          s*t*(2*Power(mrho,4) + Power(s,2) + 3*s*t + Power(t,2) - 3*Power(mrho,2)*(s + t)) - 
          2*Power(mpion,2)*(Power(mrho,4)*(s + t) + 2*s*t*(s + t) - Power(mrho,2)*(Power(s,2) + 4*s*t + Power(t,2)))
          ) - 2*eta1*eta2*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 2*Power(mpion,4)*s*t + 
          s*t*(Power(s,2) + 3*s*t + Power(t,2) - 2*Power(mrho,2)*(s + t)) + 
          Power(mpion,2)*(-4*s*t*(s + t) + Power(mrho,2)*(Power(s,2) + 4*s*t + Power(t,2))))))/
   (16.*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2))*(-Power(ma1,2) + t));

   double MaMc = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (-(eta2*(Power(mpion,2) + s)*(-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*s) + 
            s*(-Power(mrho,2) + s + 2*t))) + eta1*
        (-4*Power(mpion,6) + s*(-Power(mrho,2) + s)*(-Power(mrho,2) + s + t) + 
          Power(mpion,4)*(3*Power(mrho,2) + s + t) - 
          Power(mpion,2)*(Power(mrho,4) + 2*s*(s - t) + Power(mrho,2)*(-s + t)))))/
   (8.*(-Power(mpion,2) + s)*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

   double MaMd = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (-(eta2*(Power(mpion,2) + s)) + eta1*(-Power(mrho,2) + s + t))*
     (-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*t) + t*(-Power(mrho,2) + 2*s + t)))/
   (8.*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2))*(-Power(mpion,2) + t));

   double MaMe = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (eta2*(-2*Power(mpion,4)*(delta - 4*C4*Power(mrho,2))*(Power(mrho,2) + 4*s) + 
          Power(mpion,2)*(-2*Power(mrho,4)*(-2 + delta + 8*C4*s) + 8*delta*s*(s + t) - 
             Power(mrho,2)*((-10 + delta)*s - (-2 + delta)*t + 32*C4*s*(s + t))) + 
          s*(2*Power(mrho,4)*(-2 + delta + 4*C4*s) - 2*delta*Power(s + t,2) + 
             Power(mrho,2)*((-6 + delta)*s + (-2 + delta)*t + 8*C4*Power(s + t,2)))) + 
       eta1*(4*Power(mpion,4)*(6*C4*Power(mrho,4) + 2*delta*s + Power(mrho,2)*(1 - 2*delta - 8*C4*s)) + 
          2*delta*s*Power(s + t,2) - Power(mrho,2)*
           ((-6 + 5*delta)*Power(s,2) + 2*(-2 + 3*delta)*s*t + (-2 + delta)*Power(t,2) + 8*C4*s*Power(s + t,2)) + 
          Power(mrho,4)*((-2 + delta)*(3*s + t) + 8*C4*s*(s + 2*t)) - 
          2*Power(mpion,2)*(4*delta*s*(s + t) - 
             Power(mrho,2)*(-6*s + 7*delta*s - 2*t + 3*delta*t + 16*C4*s*(s + t)) + 
             2*Power(mrho,4)*(-2 + delta + 4*C4*(2*s + t))))))/
   (8.*Power(mrho,2)*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

   double MbMb = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) - Power(mpion,4)*(Power(mrho,4) + 2*Power(mrho,2)*t + 2*(s - 2*t)*t) + 
          Power(t,2)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 2*Power(s,2) + 2*s*t + Power(t,2)) + 
          2*Power(mpion,2)*t*(Power(mrho,4) + Power(mrho,2)*(s + t) - 2*t*(s + t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(Power(mrho,4) + 4*Power(mrho,2)*t - 2*(s - 2*t)*t) + 
          Power(t,2)*(Power(mrho,4) + 2*Power(s,2) - 2*Power(mrho,2)*(s - t) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(Power(mrho,4) - Power(mrho,2)*(s - 2*t) + 2*t*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(3*Power(mrho,4) + 2*Power(mrho,2)*(s - 3*t) - 2*(s - 2*t)*t) + 
          t*(-Power(mrho,2) + t)*(2*Power(s,2) + 2*s*t + Power(t,2) - Power(mrho,2)*(2*s + t)) - 
          2*Power(mpion,2)*(-Power(mrho,2) + t)*(2*t*(s + t) - Power(mrho,2)*(s + 2*t)))))/
   (32.*Power(-Power(ma1,2) + t,2));

   double MbMc = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*(-(eta2*(Power(mpion,2) + t)) + eta1*(-Power(mrho,2) + s + t))*
     (-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*s) + s*(-Power(mrho,2) + s + 2*t)))/
   (16.*(-Power(mpion,2) + s)*(-Power(ma1,2) + t));

   double MbMd = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
     (-(eta2*(Power(mpion,2) + t)*(-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*t) + 
            t*(-Power(mrho,2) + 2*s + t))) + eta1*
        (-4*Power(mpion,6) + t*(-Power(mrho,2) + t)*(-Power(mrho,2) + s + t) + 
          Power(mpion,4)*(3*Power(mrho,2) + s + t) + 
          Power(mpion,2)*(-Power(mrho,4) + 2*(s - t)*t + Power(mrho,2)*(-s + t)))))/
   (16.*(-Power(ma1,2) + t)*(-Power(mpion,2) + t));

   double MbMe = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(eta2*(-2*Power(mpion,4)*(delta - 4*C4*Power(mrho,2))*
           (Power(mrho,2) + 4*t) + Power(mpion,2)*
           (8*delta*t*(s + t) - 2*Power(mrho,4)*(-2 + delta + 8*C4*t) - 
             Power(mrho,2)*(-((-2 + delta)*s) + (-10 + delta)*t + 32*C4*t*(s + t))) + 
          t*(-2*delta*Power(s + t,2) + 2*Power(mrho,4)*(-2 + delta + 4*C4*t) + 
             Power(mrho,2)*((-2 + delta)*s + (-6 + delta)*t + 8*C4*Power(s + t,2)))) + 
       eta1*(2*delta*t*Power(s + t,2) - Power(mrho,2)*
           ((-2 + delta)*Power(s,2) + 2*(-2 + 3*delta)*s*t + (-6 + 5*delta)*Power(t,2) + 8*C4*t*Power(s + t,2)) + 
          Power(mrho,4)*(8*C4*t*(2*s + t) + (-2 + delta)*(s + 3*t)) + 
          4*Power(mpion,4)*(6*C4*Power(mrho,4) + 2*delta*t + Power(mrho,2)*(1 - 2*delta - 8*C4*t)) - 
          2*Power(mpion,2)*(4*delta*t*(s + t) - 
             Power(mrho,2)*(-2*s + 3*delta*s - 6*t + 7*delta*t + 16*C4*t*(s + t)) + 
             2*Power(mrho,4)*(-2 + delta + 4*C4*(s + 2*t))))))/(16.*Power(mrho,2)*(-Power(ma1,2) + t));
    
    double McMc = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + s,2) - 2*Power(mpion,2)*(Power(mrho,2) + s)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + s,2));

    double McMd = (Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*(Power(mrho,6) - 2*Power(mrho,4)*(s + t) + s*t*(s + t) + 
       Power(mpion,4)*(-Power(mrho,2) + s + t) + Power(mrho,2)*(Power(s,2) + s*t + Power(t,2)) - 
       Power(mpion,2)*(2*Power(mrho,4) - 3*Power(mrho,2)*(s + t) + Power(s + t,2))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + s)*(-Power(mpion,2) + t));

    double McMe = -(Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(mpion,4)*(-2 + 3*delta - 8*C4*Power(mrho,2)) - 2*s*t + delta*s*t + 
        2*delta*Power(t,2) + Power(mrho,4)*(-2 + delta + 8*C4*(t - u)) - 2*delta*Power(u,2) + 
        Power(mpion,2)*(8*C4*Power(mrho,4) + 2*s - delta*s + 2*t - 5*delta*t + 
           4*Power(mrho,2)*(-3 + delta + 4*C4*(t - u)) + 4*delta*u) + 
        Power(mrho,2)*(2*s - delta*s + 2*t - 3*delta*t + 4*u - 8*C4*(Power(t,2) - Power(u,2)))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + s));

    double MdMd = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + t,2) - 2*Power(mpion,2)*(Power(mrho,2) + t)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + t,2));

    double MdMe = -(Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(mpion,4)*(-2 + 3*delta - 8*C4*Power(mrho,2)) + 2*delta*Power(s,2) - 
        2*s*t + delta*s*t + Power(mrho,4)*(-2 + delta + 8*C4*(s - u)) - 2*delta*Power(u,2) + 
        Power(mpion,2)*(8*C4*Power(mrho,4) + 2*s - 5*delta*s + 2*t - delta*t + 
           4*Power(mrho,2)*(-3 + delta + 4*C4*(s - u)) + 4*delta*u) + 
        Power(mrho,2)*(2*s - 3*delta*s + 2*t - delta*t + 4*u - 8*C4*(Power(s,2) - Power(u,2)))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + t));

    double MeMe = (Power(C,2)*Power(ghat,4)*(32*Power(C4,2)*Power(mrho,8) + 2*Power(delta,2)*Power(u,2) - 
       8*C4*Power(mrho,6)*(6 - delta + 8*C4*u) - 2*delta*Power(mrho,2)*u*(6 - delta + 8*C4*u) + 
       Power(mrho,4)*(12 - Power(delta,2) + 8*C4*(6 + delta)*u + 32*Power(C4,2)*Power(u,2))))/(4.*Power(mrho,4));

    double result = MaMa + MaMb + MaMc + MaMd + MaMe + MbMb + 2*(MbMc + MbMd + MbMe) 
                    + McMc + 2*(McMd + McMe) + MdMd + 2*MdMe + MeMe;
    return(result);
}

double C1p2(double s, double t)
{
   double C = 0.059;
   double ghat = g_tilde;
  
   double u = -s - t + 2*mpion*mpion + mrho*mrho;

   double MaMa = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) + Power(mpion,4)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 4*Power(s,2) - 2*s*t) + 
          Power(s,2)*(-Power(mrho,4) + Power(s,2) - 2*Power(mrho,2)*t + 2*s*t + 2*Power(t,2)) + 
          2*Power(mpion,2)*s*(Power(mrho,4) + Power(mrho,2)*(s + t) - 2*s*(s + t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(Power(mrho,4) + 4*Power(mrho,2)*s + 4*Power(s,2) - 2*s*t) + 
          Power(s,2)*(Power(mrho,4) + Power(s,2) + 2*Power(mrho,2)*(s - t) + 2*s*t + 2*Power(t,2)) - 
          2*Power(mpion,2)*s*(Power(mrho,4) + Power(mrho,2)*(2*s - t) + 2*s*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) - 
          2*Power(mpion,2)*(-Power(mrho,2) + s)*(2*s*(s + t) - Power(mrho,2)*(2*s + t)) + 
          Power(mpion,4)*(3*Power(mrho,4) + 2*s*(2*s - t) + Power(mrho,2)*(-6*s + 2*t)) + 
          s*(-Power(mrho,2) + s)*(Power(s,2) + 2*s*t + 2*Power(t,2) - Power(mrho,2)*(s + 2*t)))))/
   (32.*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

   double MaMb = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (-(eta2*(Power(mpion,2) + s)*(-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*s) + 
            s*(-Power(mrho,2) + s + 2*t))) + eta1*
        (-4*Power(mpion,6) + s*(-Power(mrho,2) + s)*(-Power(mrho,2) + s + t) + 
          Power(mpion,4)*(3*Power(mrho,2) + s + t) - 
          Power(mpion,2)*(Power(mrho,4) + 2*s*(s - t) + Power(mrho,2)*(-s + t)))))/
   (8.*(-Power(mpion,2) + s)*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

    double MaMc = -(Power(C,2)*Power(ghat,4)*(-Power(ma1,2) + s)*
      (-2*delta*Power(mpion,2) - (-2 + delta)*Power(mrho,2) + delta*(s + t))*
      (Power(eta1,2)*(8*Power(mpion,6) + Power(s,3) - 2*Power(mrho,4)*(s - t) + Power(s,2)*t + 5*s*Power(t,2) + 
           Power(t,3) - 2*Power(mpion,2)*(s + t)*(-2*Power(mrho,2) + s + t) - 
           2*Power(mpion,4)*(2*Power(mrho,2) + s + 3*t) + Power(mrho,2)*(Power(s,2) - 2*s*t - 3*Power(t,2))) + 
        Power(eta2,2)*(4*Power(mpion,4)*(s - t) + s*(s - t)*(4*Power(mrho,2) + s - t) + 
           Power(mpion,2)*(-3*Power(s,2) + 2*s*t + Power(t,2))) + 
        eta1*eta2*(-8*Power(mpion,6) - 2*Power(s,3) - 2*Power(mpion,4)*(-2*Power(mrho,2) + s - 5*t) + 
           2*Power(mrho,4)*(s - t) + Power(s,2)*t - 6*s*Power(t,2) - Power(t,3) + 
           Power(mrho,2)*(-5*Power(s,2) + 6*s*t + 3*Power(t,2)) + 
           Power(mpion,2)*(5*Power(s,2) + 2*s*t + Power(t,2) - 4*Power(mrho,2)*(s + t)))))/
   (16.*Power(mrho,2)*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2))*
     (-2*Power(mpion,2) + s + t));

   double MaMd = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (-(eta2*(2*Power(mpion,4)*(-4*C4*Power(mrho,4) + 6*delta*s + Power(mrho,2)*(delta - 16*C4*s) - 2*delta*t) + 
            Power(mpion,2)*(2*Power(mrho,4)*(-2 + delta + 8*C4*s) + delta*(-11*Power(s,2) - 6*s*t + Power(t,2)) + 
               Power(mrho,2)*((-10 + delta)*s - (-2 + delta)*t + 32*C4*s*(s + t))) + 
            s*(-2*Power(mrho,4)*(-2 + delta + 4*C4*s) + delta*(3*Power(s,2) + 2*s*t + 3*Power(t,2)) + 
               Power(mrho,2)*(3*(2 + delta)*s + (2 - 5*delta)*t - 8*C4*Power(s + t,2))))) + 
       eta1*(8*delta*Power(mpion,6) + 2*Power(mpion,4)*
           (12*C4*Power(mrho,4) - 2*Power(mrho,2)*(-1 + 3*delta + 8*C4*s) + 3*delta*(s - t)) + 
          delta*(3*Power(s,3) + 5*Power(s,2)*t + 7*s*Power(t,2) + Power(t,3)) - 
          2*Power(mrho,2)*((-3 + 2*delta)*Power(s,2) + 2*(-1 + 2*delta)*s*t + (-1 + 2*delta)*Power(t,2) + 
             4*C4*s*Power(s + t,2)) + Power(mrho,4)*((-6 + delta)*s + (-2 + 3*delta)*t + 8*C4*s*(s + 2*t)) - 
          2*Power(mpion,2)*(delta*(5*Power(s,2) + 6*s*t + Power(t,2)) - 
             Power(mrho,2)*(-6*s + 9*delta*s - 2*t + 5*delta*t + 16*C4*s*(s + t)) + 
             2*Power(mrho,4)*(-2 + delta + 4*C4*(2*s + t))))))/
   (16.*Power(mrho,2)*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

   double MbMb = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + s,2) - 2*Power(mpion,2)*(Power(mrho,2) + s)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + s,2));

   double MbMc = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(-2*Power(mrho,2) + delta*u)*
     (-4*Power(mpion,4)*(-5*Power(mrho,2) + u) + u*(-3*Power(mrho,4) - s*u + Power(mrho,2)*(3*s + 2*t + u)) + 
       Power(mpion,2)*(12*Power(mrho,4) + u*(4*s + u) - Power(mrho,2)*(8*s + 12*t + 9*u))))/
   (16.*Power(mrho,4)*(-Power(mpion,2) + s)*(-Power(mrho,2) + u));

   double MbMd = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(6*delta*Power(mpion,6) + Power(mrho,6)*(-2 + 3*delta + 8*C4*s) + 
       delta*s*t*(s + t) + Power(mrho,2)*(3*delta*Power(s,2) + (2 + 3*delta)*s*t + delta*Power(t,2)) - 
       Power(mpion,4)*((-2 + 9*delta)*Power(mrho,2) - 8*C4*Power(mrho,4) + delta*(9*s + t)) - 
       2*Power(mrho,4)*(-s + 3*delta*s - t + 2*delta*t + 4*C4*s*(s + 2*t)) + 
       Power(mpion,2)*(-8*C4*Power(mrho,6) - 2*Power(mrho,4)*(-2 + delta - 8*C4*s) - 
          Power(mrho,2)*(2*s + 5*delta*s + 2*t - 7*delta*t) + delta*(3*Power(s,2) - Power(t,2)))))/
   (16.*Power(mrho,4)*(-Power(mpion,2) + s));

   double McMc = (Power(C,2)*Power(ghat,4)*Power(-2*Power(mrho,2) + delta*u,2)*
     (-3*Power(mrho,4)*u + Power(u,3) + 2*Power(mrho,2)*(Power(s,2) - 2*s*t + Power(t,2) - Power(u,2)) - 
       4*Power(mpion,2)*(-3*Power(mrho,4) - 2*Power(mrho,2)*u + Power(u,2))))/
   (16.*Power(mrho,6)*Power(-Power(mrho,2) + u,2));

   double McMd = (Power(C,2)*Power(ghat,4)*(-2*Power(mrho,2) + delta*u)*
     (24*C4*Power(mrho,6)*(s - t) + Power(mrho,4)*
        ((-10 + delta)*s + 10*t - delta*t + 6*delta*u + 8*C4*(-s + t)*u) + 
       delta*u*(Power(s,2) - Power(t,2) - 2*Power(u,2)) - 
       2*delta*Power(mpion,2)*(12*Power(mrho,4) + (s - t - 4*u)*u + Power(mrho,2)*(-s + t + 8*u)) + 
       Power(mrho,2)*(-5*delta*Power(s,2) - 3*delta*Power(t,2) + (-2 + delta)*t*u + 4*delta*Power(u,2) + 
          s*(8*delta*t - (-2 + delta)*u))))/(32.*Power(mrho,6)*(-Power(mrho,2) + u));

   double MdMd = (Power(C,2)*Power(ghat,4)*(32*Power(C4,2)*Power(mrho,10) + 4*Power(delta,2)*Power(mpion,4)*u + 
       2*Power(delta,2)*t*(s + t)*u - 8*C4*Power(mrho,8)*(6 - 3*delta + 8*C4*u) - 
       2*delta*Power(mrho,4)*((-8 + delta)*s + 2*t + 2*u + delta*u - 4*C4*(3*s + t)*u) - 
       2*delta*Power(mpion,2)*(-16*C4*Power(mrho,6) + delta*u*(s + 3*t + 2*u) + 
          Power(mrho,4)*(6 - 7*delta + 16*C4*u) + Power(mrho,2)*(delta*(s - t) + (2 - 7*delta)*u)) + 
       Power(mrho,6)*(12 - 8*delta + Power(delta,2) + 32*Power(C4,2)*Power(u,2) - 
          8*C4*(delta*(5*s - t) + 3*(-2 + delta)*u)) + 
       delta*Power(mrho,2)*(delta*Power(t,2) + delta*s*(3*s - 2*u) - 4*t*(delta*s + (-1 + delta)*u))))/
   (16.*Power(mrho,6));

    double result = MaMa + MaMb + MaMc + MaMd + MbMb + 2*(MbMc + MbMd) + McMc + 2*McMd + MdMd;
    return(result);
}

double C1p3(double s, double t)
{
   double C = 0.059;
   double ghat = g_tilde;
  
   double u = -s - t + 2*mpion*mpion + mrho*mrho;

   double MaMa = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) - Power(mpion,4)*(Power(mrho,4) + 2*Power(mrho,2)*t + 2*(s - 2*t)*t) + 
          Power(t,2)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 2*Power(s,2) + 2*s*t + Power(t,2)) + 
          2*Power(mpion,2)*t*(Power(mrho,4) + Power(mrho,2)*(s + t) - 2*t*(s + t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(Power(mrho,4) + 4*Power(mrho,2)*t - 2*(s - 2*t)*t) + 
          Power(t,2)*(Power(mrho,4) + 2*Power(s,2) - 2*Power(mrho,2)*(s - t) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(Power(mrho,4) - Power(mrho,2)*(s - 2*t) + 2*t*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(3*Power(mrho,4) + 2*Power(mrho,2)*(s - 3*t) - 2*(s - 2*t)*t) + 
          t*(-Power(mrho,2) + t)*(2*Power(s,2) + 2*s*t + Power(t,2) - Power(mrho,2)*(2*s + t)) - 
          2*Power(mpion,2)*(-Power(mrho,2) + t)*(2*t*(s + t) - Power(mrho,2)*(s + 2*t)))))/
   (32.*Power(-Power(ma1,2) + t,2));

   double MaMb = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
     (-(eta2*(Power(mpion,2) + t)*(-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*t) + 
            t*(-Power(mrho,2) + 2*s + t))) + eta1*
        (-4*Power(mpion,6) + t*(-Power(mrho,2) + t)*(-Power(mrho,2) + s + t) + 
          Power(mpion,4)*(3*Power(mrho,2) + s + t) + 
          Power(mpion,2)*(-Power(mrho,4) + 2*(s - t)*t + Power(mrho,2)*(-s + t)))))/
   (8.*(-Power(ma1,2) + t)*(-Power(mpion,2) + t));

   double MaMc = -(Power(C,2)*Power(ghat,4)*(-2*delta*Power(mpion,2) - (-2 + delta)*Power(mrho,2) + delta*(s + t))*
      (Power(eta2,2)*(-4*Power(mpion,4)*(s - t) + (s - t)*(-4*Power(mrho,2) + s - t)*t + 
           Power(mpion,2)*(Power(s,2) + 2*s*t - 3*Power(t,2))) + 
        Power(eta1,2)*(8*Power(mpion,6) + Power(s,3) + 2*Power(mrho,4)*(s - t) + 5*Power(s,2)*t + s*Power(t,2) + 
           Power(t,3) - 2*Power(mpion,2)*(s + t)*(-2*Power(mrho,2) + s + t) - 
           2*Power(mpion,4)*(2*Power(mrho,2) + 3*s + t) + Power(mrho,2)*(-3*Power(s,2) - 2*s*t + Power(t,2))) + 
        eta1*eta2*(-8*Power(mpion,6) - Power(s,3) - 2*Power(mrho,4)*(s - t) + 
           2*Power(mpion,4)*(2*Power(mrho,2) + 5*s - t) - 6*Power(s,2)*t + s*Power(t,2) - 2*Power(t,3) + 
           Power(mrho,2)*(3*Power(s,2) + 6*s*t - 5*Power(t,2)) + 
           Power(mpion,2)*(Power(s,2) + 2*s*t + 5*Power(t,2) - 4*Power(mrho,2)*(s + t)))))/
   (16.*Power(mrho,2)*(-Power(ma1,2) + t)*(-2*Power(mpion,2) + s + t));

   double MaMd = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(-(eta2*
          (-2*Power(mpion,4)*(4*C4*Power(mrho,4) + 2*delta*(s - 3*t) - Power(mrho,2)*(delta - 16*C4*t)) + 
            Power(mpion,2)*(2*Power(mrho,4)*(-2 + delta + 8*C4*t) + delta*(Power(s,2) - 6*s*t - 11*Power(t,2)) + 
               Power(mrho,2)*(-((-2 + delta)*s) + (-10 + delta)*t + 32*C4*t*(s + t))) + 
            t*(-2*Power(mrho,4)*(-2 + delta + 4*C4*t) + delta*(3*Power(s,2) + 2*s*t + 3*Power(t,2)) + 
               Power(mrho,2)*((2 - 5*delta)*s + 3*(2 + delta)*t - 8*C4*Power(s + t,2))))) + 
       eta1*(8*delta*Power(mpion,6) + delta*(Power(s,3) + 7*Power(s,2)*t + 5*s*Power(t,2) + 3*Power(t,3)) - 
          2*Power(mrho,2)*((-1 + 2*delta)*Power(s,2) + 2*(-1 + 2*delta)*s*t + (-3 + 2*delta)*Power(t,2) + 
             4*C4*t*Power(s + t,2)) + Power(mrho,4)*((-2 + 3*delta)*s + (-6 + delta)*t + 8*C4*t*(2*s + t)) + 
          Power(mpion,4)*(24*C4*Power(mrho,4) + 6*delta*(-s + t) - 4*Power(mrho,2)*(-1 + 3*delta + 8*C4*t)) - 
          2*Power(mpion,2)*(delta*(Power(s,2) + 6*s*t + 5*Power(t,2)) - 
             Power(mrho,2)*(-2*s + 5*delta*s - 6*t + 9*delta*t + 16*C4*t*(s + t)) + 
             2*Power(mrho,4)*(-2 + delta + 4*C4*(s + 2*t))))))/(16.*Power(mrho,2)*(-Power(ma1,2) + t));

   double MbMb = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + t,2) - 2*Power(mpion,2)*(Power(mrho,2) + t)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + t,2));

   double MbMc = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(-2*Power(mrho,2) + delta*u)*
     (-4*Power(mpion,4)*(-5*Power(mrho,2) + u) + u*(-3*Power(mrho,4) - t*u + Power(mrho,2)*(2*s + 3*t + u)) + 
       Power(mpion,2)*(12*Power(mrho,4) + u*(4*t + u) - Power(mrho,2)*(12*s + 8*t + 9*u))))/
   (16.*Power(mrho,4)*(-Power(mpion,2) + t)*(-Power(mrho,2) + u));

   double MbMd = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(6*delta*Power(mpion,6) + delta*s*t*(s + t) + 
       Power(mrho,6)*(-2 + 3*delta + 8*C4*t) + 
       Power(mrho,2)*(delta*Power(s,2) + (2 + 3*delta)*s*t + 3*delta*Power(t,2)) - 
       2*Power(mrho,4)*(-s + 2*delta*s - t + 3*delta*t + 4*C4*t*(2*s + t)) - 
       Power(mpion,4)*((-2 + 9*delta)*Power(mrho,2) - 8*C4*Power(mrho,4) + delta*(s + 9*t)) - 
       Power(mpion,2)*(8*C4*Power(mrho,6) + 2*Power(mrho,4)*(-2 + delta - 8*C4*t) + 
          Power(mrho,2)*(2*s - 7*delta*s + 2*t + 5*delta*t) + delta*(Power(s,2) - 3*Power(t,2)))))/
   (16.*Power(mrho,4)*(-Power(mpion,2) + t));

   double McMc = (Power(C,2)*Power(ghat,4)*Power(-2*Power(mrho,2) + delta*u,2)*
     (-3*Power(mrho,4)*u + Power(u,3) + 2*Power(mrho,2)*(Power(s,2) - 2*s*t + Power(t,2) - Power(u,2)) - 
       4*Power(mpion,2)*(-3*Power(mrho,4) - 2*Power(mrho,2)*u + Power(u,2))))/
   (16.*Power(mrho,6)*Power(-Power(mrho,2) + u,2));

   double McMd = -(Power(C,2)*Power(ghat,4)*(-2*Power(mrho,2) + delta*u)*
      (24*C4*Power(mrho,6)*(s - t) + Power(mrho,4)*
         ((-10 + delta)*s + 10*t - delta*t - 6*delta*u + 8*C4*(-s + t)*u) + 
        delta*u*(Power(s,2) - Power(t,2) + 2*Power(u,2)) - 
        2*delta*Power(mpion,2)*(-12*Power(mrho,4) + Power(mrho,2)*(-s + t - 8*u) + u*(s - t + 4*u)) + 
        Power(mrho,2)*(3*delta*Power(s,2) + 5*delta*Power(t,2) + (-2 + delta)*t*u - 4*delta*Power(u,2) - 
           s*(8*delta*t + (-2 + delta)*u))))/(32.*Power(mrho,6)*(-Power(mrho,2) + u));

   double MdMd = (Power(C,2)*Power(ghat,4)*(32*Power(C4,2)*Power(mrho,10) + 4*Power(delta,2)*Power(mpion,4)*u + 
       2*Power(delta,2)*s*(s + t)*u - 8*C4*Power(mrho,8)*(6 - 3*delta + 8*C4*u) - 
       2*delta*Power(mrho,4)*(2*s - 8*t + delta*t + 2*u + delta*u - 4*C4*(s + 3*t)*u) + 
       delta*Power(mrho,2)*(delta*Power(s,2) + delta*t*(3*t - 2*u) - 4*s*(delta*t + (-1 + delta)*u)) - 
       2*delta*Power(mpion,2)*(-16*C4*Power(mrho,6) + delta*u*(3*s + t + 2*u) + 
          Power(mrho,4)*(6 - 7*delta + 16*C4*u) + Power(mrho,2)*(-(delta*s) + delta*t + 2*u - 7*delta*u)) + 
       Power(mrho,6)*(12 - 8*delta + Power(delta,2) + 32*Power(C4,2)*Power(u,2) - 
          8*C4*(-(delta*s) + 5*delta*t - 6*u + 3*delta*u))))/(16.*Power(mrho,6));
   
   double result = MaMa + (MaMb + MaMc + MaMd) + MbMb + 2*(MbMc + MbMd) + McMc + 2*McMd + MdMd;

   return(result);
}

/*************************************************************************************/
//pi + rho -> omega -> pi + gamma  (C1.4 + C1.5 + C1.6)
/*************************************************************************************/

double C1p4(double s, double t)
{
   double C = 0.059;
   double gorp = 22.6;
  
   double MaMa = (Power(C,2)*Power(gorp,4)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
       Power(mpion,4)*(Power(mrho,4) + 4*Power(s,2) - 2*s*t) + 
       Power(s,2)*(Power(mrho,4) + Power(s,2) + 2*s*t + 2*Power(t,2) - 2*Power(mrho,2)*(s + t)) - 
       2*Power(mpion,2)*s*(Power(mrho,4) + 2*s*(s + t) - Power(mrho,2)*(2*s + t))))/
   (8.*Power(-Power(momega,2) + s,2));

   double MaMb = (Power(C,2)*Power(gorp,4)*(Power(mpion,8) - 4*Power(mpion,2)*s*t*(-Power(mrho,2) + s + t) + 
       Power(mpion,4)*(2*s*t - Power(mrho,2)*(s + t)) + 
       s*t*(Power(s,2) + 3*s*t + Power(t,2) - Power(mrho,2)*(s + t))))/
   (8.*(-Power(momega,2) + s)*(-Power(momega,2) + t));

   double MbMb = (Power(C,2)*Power(gorp,4)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
       Power(mpion,4)*(Power(mrho,4) - 2*(s - 2*t)*t) + 
       Power(t,2)*(Power(mrho,4) + 2*Power(s,2) + 2*s*t + Power(t,2) - 2*Power(mrho,2)*(s + t)) - 
       2*Power(mpion,2)*t*(Power(mrho,4) + 2*t*(s + t) - Power(mrho,2)*(s + 2*t))))/
   (8.*Power(-Power(momega,2) + t,2));

   double result = MaMa + 2*MaMb + MbMb;
   return(result);
}

double C1p5(double s, double t)
{
   double C = 0.059;
   double gorp = 22.6;
  
   double MaMa = (Power(C,2)*Power(gorp,4)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
       Power(mpion,4)*(Power(mrho,4) + 4*Power(s,2) - 2*s*t) + 
       Power(s,2)*(Power(mrho,4) + Power(s,2) + 2*s*t + 2*Power(t,2) - 2*Power(mrho,2)*(s + t)) - 
       2*Power(mpion,2)*s*(Power(mrho,4) + 2*s*(s + t) - Power(mrho,2)*(2*s + t))))/
   (8.*Power(-Power(momega,2) + s,2));

   double result = MaMa;
   return(result);
}

double C1p6(double s, double t)
{
   double C = 0.059;
   double gorp = 22.6;
  
   double MaMa = (Power(C,2)*Power(gorp,4)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
       Power(mpion,4)*(Power(mrho,4) - 2*(s - 2*t)*t) + 
       Power(t,2)*(Power(mrho,4) + 2*Power(s,2) + 2*s*t + Power(t,2) - 2*Power(mrho,2)*(s + t)) - 
       2*Power(mpion,2)*t*(Power(mrho,4) + 2*t*(s + t) - Power(mrho,2)*(s + 2*t))))/
   (8.*Power(-Power(momega,2) + t,2));

   double result = MaMa;
   return(result);
}

/*************************************************************************************/
//pi + pi -> rho + gamma  (C2.1 + C2.2)
/*************************************************************************************/

double C2p1(double s, double t)
{
   double C = 0.059;
   double ghat = g_tilde;
	
   double u = -s - t + 2*mpion*mpion + mrho*mrho; 

   double MaMa = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (Power(eta1,2)*(Power(mpion,8) + Power(mpion,6)*(2*Power(mrho,2) - 4*t) + Power(mrho,4)*Power(s + t,2) + 
          Power(s + t,2)*(Power(s,2) + Power(t,2)) - 
          2*Power(mrho,2)*(Power(s,3) + 2*Power(s,2)*t + 2*s*Power(t,2) + Power(t,3)) + 
          Power(mpion,4)*(5*Power(mrho,4) + 4*Power(s,2) + 2*s*t + 6*Power(t,2) - 2*Power(mrho,2)*(5*s + 3*t)) - 
          2*Power(mpion,2)*(2*Power(mrho,4)*(s + t) - Power(mrho,2)*(4*Power(s,2) + 5*s*t + 3*Power(t,2)) + 
             2*(Power(s,3) + Power(s,2)*t + s*Power(t,2) + Power(t,3)))) - 
       2*eta1*eta2*(Power(mpion,8) - 4*Power(mpion,6)*(-Power(mrho,2) + t) + 
          Power(-Power(mrho,2) + s + t,2)*(Power(s,2) + Power(t,2) - 2*Power(mrho,2)*(s + t)) + 
          Power(mpion,4)*(9*Power(mrho,4) + 4*Power(s,2) + 2*s*t + 6*Power(t,2) - 2*Power(mrho,2)*(7*s + 6*t)) - 
          2*Power(mpion,2)*(-2*Power(mrho,6) + Power(mrho,4)*(7*s + 6*t) - 
             Power(mrho,2)*(7*Power(s,2) + 9*s*t + 6*Power(t,2)) + 
             2*(Power(s,3) + Power(s,2)*t + s*Power(t,2) + Power(t,3)))) + 
       Power(eta2,2)*(Power(mpion,8) + Power(mpion,6)*(6*Power(mrho,2) - 4*t) + 
          Power(-Power(mrho,2) + s + t,2)*(4*Power(mrho,4) + Power(s,2) + Power(t,2) - 4*Power(mrho,2)*(s + t)) + 
          Power(mpion,4)*(17*Power(mrho,4) + 4*Power(s,2) + 2*s*t + 6*Power(t,2) - 2*Power(mrho,2)*(10*s + 9*t)) - 
          2*Power(mpion,2)*(-7*Power(mrho,6) + Power(mrho,4)*(15*s + 14*t) - 
             Power(mrho,2)*(10*Power(s,2) + 15*s*t + 9*Power(t,2)) + 
             2*(Power(s,3) + Power(s,2)*t + s*Power(t,2) + Power(t,3))))))/
   (32.*(Power(ma1,4) + Power(-2*Power(mpion,2) - Power(mrho,2) + s + t,2) + 
       Power(ma1,2)*(Power(Gammaa1,2) - 4*Power(mpion,2) - 2*Power(mrho,2) + 2*s + 2*t)));

   double MaMb = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*(Power(ma1,2) - 2*Power(mpion,2) - Power(mrho,2) + s + t)*
     (Power(eta2,2)*(-Power(mpion,8) + Power(mpion,6)*(-2*Power(mrho,2) + 4*t) - 
          2*Power(mpion,2)*t*(Power(mrho,4) + Power(s,2) - 2*s*t - 2*Power(t,2) + Power(mrho,2)*(-2*s + 3*t)) - 
          Power(mpion,4)*(Power(mrho,4) + 2*t*(s + 3*t) - Power(mrho,2)*(s + 6*t)) + 
          t*(-2*Power(mrho,6) + Power(s,3) - 2*s*Power(t,2) - Power(t,3) + Power(mrho,4)*(5*s + t) - 
             Power(mrho,2)*(4*Power(s,2) + s*t - 2*Power(t,2)))) + 
       Power(eta1,2)*(-Power(mpion,8) + Power(mpion,6)*(-2*Power(mrho,2) + 4*t) - 
          t*(-Power(s,3) + 2*s*Power(t,2) + Power(t,3) + Power(mrho,4)*(s + t) - Power(mrho,2)*t*(3*s + 2*t)) + 
          Power(mpion,4)*(-3*Power(mrho,4) - 2*t*(s + 3*t) + Power(mrho,2)*(5*s + 6*t)) - 
          2*Power(mpion,2)*(-(Power(mrho,4)*(s + t)) + t*(Power(s,2) - 2*s*t - 2*Power(t,2)) + 
             Power(mrho,2)*(Power(s,2) + 2*s*t + 3*Power(t,2)))) + 
       2*eta1*eta2*(Power(mpion,8) + Power(mpion,6)*(2*Power(mrho,2) - 4*t) + 
          2*Power(mpion,4)*(2*Power(mrho,4) + t*(s + 3*t) - Power(mrho,2)*(2*s + 3*t)) + 
          t*(-Power(mrho,6) - Power(s,3) + 2*s*Power(t,2) + Power(t,3) + Power(mrho,4)*(s + 2*t) + 
             Power(mrho,2)*(Power(s,2) - 2*s*t - 2*Power(t,2))) + 
          Power(mpion,2)*(Power(mrho,6) - 2*Power(mrho,4)*(s + 2*t) + 2*t*(Power(s,2) - 2*s*t - 2*Power(t,2)) + 
             Power(mrho,2)*(Power(s,2) + 2*s*t + 6*Power(t,2))))))/
   (16.*(-Power(ma1,2) + t)*(Power(ma1,4) + Power(-2*Power(mpion,2) - Power(mrho,2) + s + t,2) + 
       Power(ma1,2)*(Power(Gammaa1,2) - 4*Power(mpion,2) - 2*Power(mrho,2) + 2*s + 2*t)));

   double MaMc = -(Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(ma1,2) - 2*Power(mpion,2) - Power(mrho,2) + s + t)*
      (Power(eta1,2)*(2*Power(mpion,6) + Power(mpion,4)*(-2*Power(mrho,2) + 5*s - 4*t) + 
           s*(s + t)*(-Power(mrho,2) + s + t) + 
           Power(mpion,2)*(2*Power(mrho,4) - 4*Power(s,2) + Power(mrho,2)*(s - 2*t) - 2*s*t + 2*Power(t,2))) + 
        eta1*eta2*(-5*Power(mpion,6) + Power(mrho,4)*(-s + t) - (2*s - t)*Power(s + t,2) + 
           Power(mpion,4)*(4*Power(mrho,2) - 10*s + 11*t) + Power(mrho,2)*(3*Power(s,2) + s*t - 2*Power(t,2)) + 
           Power(mpion,2)*(-Power(mrho,4) + 9*Power(s,2) + 2*s*t - 7*Power(t,2) + Power(mrho,2)*(-7*s + 6*t))) + 
        Power(eta2,2)*(3*Power(mpion,6) + Power(mpion,4)*(-2*Power(mrho,2) + 5*s - 7*t) + 
           (s - t)*Power(-Power(mrho,2) + s + t,2) - 
           Power(mpion,2)*(Power(mrho,4) + Power(mrho,2)*(-6*s + 4*t) + 5*(Power(s,2) - Power(t,2))))))/
   (8.*(-Power(mpion,2) - Power(mrho,2) + s + t)*
     (Power(ma1,4) + Power(-2*Power(mpion,2) - Power(mrho,2) + s + t,2) + 
       Power(ma1,2)*(Power(Gammaa1,2) - 4*Power(mpion,2) - 2*Power(mrho,2) + 2*s + 2*t)));

   double MaMd = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*(Power(ma1,2) - 2*Power(mpion,2) - Power(mrho,2) + s + t)*
     (-(eta1*(-2*Power(mpion,2) + s)) + eta2*(-3*Power(mpion,2) - Power(mrho,2) + s + t))*
     (Power(mpion,4) + t*(-Power(mrho,2) + 2*s + t) - Power(mpion,2)*(Power(mrho,2) + 2*t)))/
   (8.*(-Power(mpion,2) + t)*(Power(ma1,4) + Power(-2*Power(mpion,2) - Power(mrho,2) + s + t,2) + 
       Power(ma1,2)*(Power(Gammaa1,2) - 4*Power(mpion,2) - 2*Power(mrho,2) + 2*s + 2*t)));

   double MaMe = -(Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(Power(ma1,2) - 2*Power(mpion,2) - Power(mrho,2) + s + t)*
      (eta1*(Power(mpion,4)*(4*Power(mrho,2) - 8*C4*Power(mrho,4)) + 8*C4*Power(mrho,6)*(s + t) - 
           2*delta*Power(s,2)*(s + t) + Power(mrho,2)*
            ((6 + delta)*Power(s,2) + 8*s*t + 4*Power(t,2) + 8*C4*Power(s,2)*(s + t)) - 
           Power(mrho,4)*(-((-6 + delta)*s) + 4*t + 8*C4*(2*Power(s,2) + 2*s*t + Power(t,2))) + 
           2*Power(mpion,2)*(-8*C4*Power(mrho,6) + 2*delta*Power(s,2) - 
              Power(mrho,2)*((6 + delta)*s + 8*C4*Power(s,2) + 4*t) + 4*Power(mrho,4)*(1 + 2*C4*(2*s + t)))) + 
        eta2*(Power(mpion,4)*(-4*Power(mrho,2) + 8*C4*Power(mrho,4)) + 
           Power(mpion,2)*(32*C4*Power(mrho,6) - 4*delta*Power(s,2) + 
              Power(mrho,2)*(14*s + 5*delta*s + 16*C4*Power(s,2) + 8*t) + 
              Power(mrho,4)*(-18 + delta - 16*C4*(3*s + t))) - 
           (-Power(mrho,2) + s + t)*(16*C4*Power(mrho,6) - 2*delta*Power(s,2) + 
              Power(mrho,2)*(3*(2 + delta)*s + 8*C4*Power(s,2) + 4*t) + Power(mrho,4)*(-10 + delta - 8*C4*(3*s + t))
              ))))/(8.*Power(mrho,2)*(Power(ma1,4) + Power(-2*Power(mpion,2) - Power(mrho,2) + s + t,2) + 
       Power(ma1,2)*(Power(Gammaa1,2) - 4*Power(mpion,2) - 2*Power(mrho,2) + 2*s + 2*t)));

   double MbMb = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) - 4*Power(mpion,6)*t + 
          Power(t,2)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 2*Power(s,2) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(-2*Power(mrho,4) + Power(mrho,2)*s + 2*t*(s + t)) + 
          Power(mpion,4)*(-Power(mrho,4) + 2*t*(s + 3*t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*(Power(mrho,2) + 2*t) + 
          Power(t,2)*(Power(mrho,4) + 2*Power(s,2) - 2*Power(mrho,2)*(s - t) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(2*t*(s + t) + Power(mrho,2)*(s + 3*t)) + 
          Power(mpion,4)*(Power(mrho,4) + 6*Power(mrho,2)*t + 2*t*(s + 3*t))) + 
       Power(eta1,2)*(Power(mpion,8) + Power(mpion,6)*(2*Power(mrho,2) - 4*t) - 
          2*Power(mpion,2)*(-Power(mrho,2) + s + t)*(-Power(mrho,4) - Power(mrho,2)*t + 2*Power(t,2)) + 
          t*(-Power(mrho,2) + t)*(2*Power(s,2) + 2*s*t + Power(t,2) - Power(mrho,2)*(2*s + t)) + 
          Power(mpion,4)*(Power(mrho,4) - 2*Power(mrho,2)*(s + 3*t) + 2*t*(s + 3*t)))))/
   (32.*Power(-Power(ma1,2) + t,2));

   double MbMc = -(Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
      (-Power(mpion,4) + Power(s,2) - Power(t,2) + Power(mrho,2)*(-s + t) + 
        Power(mpion,2)*(Power(mrho,2) - 2*s + 2*t))*
      (eta1*(-2*Power(mpion,2) + s) - eta2*(-3*Power(mpion,2) - Power(mrho,2) + s + u)))/
   (16.*(-Power(ma1,2) + t)*(-Power(mpion,2) + u));

   double MbMd = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
     (-(eta1*(-6*Power(mpion,6) + s*t*(-Power(mrho,2) + t) + Power(mpion,4)*(-10*Power(mrho,2) + 5*s + 8*t + 4*u) + 
            Power(mpion,2)*(-2*t*(s + t + 2*u) + Power(mrho,2)*(s + 4*t + 2*u)))) + 
       eta2*(-5*Power(mpion,6) + t*(-Power(mrho,4) + t*(s - u) + Power(mrho,2)*(s + t + u)) + 
          Power(mpion,4)*(2*Power(mrho,2) + 3*(s + 2*t + u)) - 
          Power(mpion,2)*(-3*Power(mrho,4) + t*(t + 2*u) + Power(mrho,2)*(3*s + 5*t + 3*u)))))/
   (16.*(-Power(ma1,2) + t)*(-Power(mpion,2) + t));

   double MbMe = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(eta2*(8*C4*Power(mrho,6)*t - 2*delta*Power(s,2)*t + 
          Power(mrho,2)*(-4*Power(mpion,4) + 8*C4*Power(s,2)*t + (2*s + 3*delta*s - 4*t)*t + 
             Power(mpion,2)*(2*s - delta*s + 8*t)) + 
          Power(mrho,4)*(-((-2 + delta)*Power(mpion,2)) + (-6 + delta)*t - 
             8*C4*(-Power(mpion,4) + 2*Power(mpion,2)*t + (2*s - t)*t))) + 
       eta1*(2*delta*Power(s,2)*t + 8*C4*Power(mrho,6)*(-2*Power(mpion,2) + t) - 
          Power(mrho,2)*(-4*Power(mpion,4) - 2*Power(s,2) + delta*Power(s,2) + 8*C4*Power(s,2)*t - 4*Power(t,2) + 
             2*Power(mpion,2)*((2 + delta)*s + 4*t)) + 
          Power(mrho,4)*(8*Power(mpion,2) - 2*s + delta*s - 4*t - 
             8*C4*(Power(mpion,4) + Power(t,2) - 2*Power(mpion,2)*(s + t))))))/
   (16.*Power(mrho,2)*(-Power(ma1,2) + t));

   double McMc = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(s + t,2) - 2*Power(mpion,2)*(2*Power(mrho,2) + s + t)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) - Power(mrho,2) + s + t,2));

   double McMd = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*
      (-2*Power(mpion,6) + Power(mpion,4)*(6*Power(mrho,2) + 3*s + 4*t) - 
        Power(mpion,2)*(Power(s,2) + Power(mrho,2)*(5*s - 2*t) + 4*s*t + 2*Power(t,2)) + 
        s*(Power(mrho,2)*(s - t) + t*(s + t))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + t)*(-Power(mpion,2) - Power(mrho,2) + s + t));

   double McMe = -(Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(mpion,4)*(-2 + 3*delta - 8*C4*Power(mrho,2)) - 2*delta*Power(s,2) + 
        Power(mrho,4)*(-2 + delta - 8*C4*(s - t)) + 2*delta*Power(t,2) - 2*t*u + delta*t*u + 
        Power(mpion,2)*(8*C4*Power(mrho,4) + 4*delta*s + 4*Power(mrho,2)*(-3 + delta - 4*C4*(s - t)) + 2*t - 
           5*delta*t + 2*u - delta*u) + Power(mrho,2)*
         (4*s + 2*t - 3*delta*t + 8*C4*(Power(s,2) - Power(t,2)) + 2*u - delta*u)))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + u));

   double MdMd = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + t,2) - 2*Power(mpion,2)*(Power(mrho,2) + t)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + t,2));

   double MdMe = -(Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(mpion,4)*(-2 + 3*delta - 8*C4*Power(mrho,2)) - 2*delta*Power(s,2) + 
        Power(mrho,4)*(-2 + delta - 8*C4*(s - u)) - 2*t*u + delta*t*u + 2*delta*Power(u,2) + 
        Power(mpion,2)*(8*C4*Power(mrho,4) + 4*delta*s + 2*t - delta*t + 
           4*Power(mrho,2)*(-3 + delta - 4*C4*(s - u)) + 2*u - 5*delta*u) + 
        Power(mrho,2)*(4*s + 2*t - delta*t + 2*u - 3*delta*u + 8*C4*(Power(s,2) - Power(u,2)))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + t));

   double MeMe = (Power(C,2)*Power(ghat,4)*(32*Power(C4,2)*Power(mrho,8) + 2*Power(delta,2)*Power(s,2) - 
       8*C4*Power(mrho,6)*(6 - delta + 8*C4*s) - 2*delta*Power(mrho,2)*s*(6 - delta + 8*C4*s) + 
       Power(mrho,4)*(12 - Power(delta,2) + 8*C4*(6 + delta)*s + 32*Power(C4,2)*Power(s,2))))/(4.*Power(mrho,4));

   double result = MaMa + MaMb + MaMc + MaMd + MaMe + MbMb + 2*(MbMc + MbMd + MbMe) 
                    + McMc + 2*(McMd + McMe) + MdMd + 2*MdMe + MeMe;
   return(result);
}

double C2p2(double s, double t)
{
   double C = 0.059;
   double ghat = g_tilde;
	
   double u = -s - t + 2*mpion*mpion + mrho*mrho; 
   double MaMa = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) - 4*Power(mpion,6)*t + 
          Power(t,2)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 2*Power(s,2) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(-2*Power(mrho,4) + Power(mrho,2)*s + 2*t*(s + t)) + 
          Power(mpion,4)*(-Power(mrho,4) + 2*t*(s + 3*t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*(Power(mrho,2) + 2*t) + 
          Power(t,2)*(Power(mrho,4) + 2*Power(s,2) - 2*Power(mrho,2)*(s - t) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(2*t*(s + t) + Power(mrho,2)*(s + 3*t)) + 
          Power(mpion,4)*(Power(mrho,4) + 6*Power(mrho,2)*t + 2*t*(s + 3*t))) + 
       Power(eta1,2)*(Power(mpion,8) + Power(mpion,6)*(2*Power(mrho,2) - 4*t) - 
          2*Power(mpion,2)*(-Power(mrho,2) + s + t)*(-Power(mrho,4) - Power(mrho,2)*t + 2*Power(t,2)) + 
          t*(-Power(mrho,2) + t)*(2*Power(s,2) + 2*s*t + Power(t,2) - Power(mrho,2)*(2*s + t)) + 
          Power(mpion,4)*(Power(mrho,4) - 2*Power(mrho,2)*(s + 3*t) + 2*t*(s + 3*t)))))/
   (32.*Power(-Power(ma1,2) + t,2));

   double MaMb = -(Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
      (-(eta2*(Power(mpion,2) + t)*(Power(mpion,4) + t*(-Power(mrho,2) + 2*s + t) - 
             Power(mpion,2)*(Power(mrho,2) + 2*t))) + 
        eta1*(2*Power(mpion,6) + Power(mpion,4)*(-2*Power(mrho,2) + s - 4*t) + s*t*(-Power(mrho,2) + t) + 
           Power(mpion,2)*(2*Power(mrho,4) + 2*t*(s + t) - Power(mrho,2)*(s + 2*t)))))/
   (8.*(-Power(ma1,2) + t)*(-Power(mpion,2) + t));

   double MaMc = (Power(C,2)*Power(ghat,4)*(2*Power(mrho,2) - delta*s)*
     (eta1*eta2*(8*Power(mpion,6) + Power(s,3) - 3*Power(s,2)*t - 10*s*Power(t,2) - 8*Power(t,3) + 
          Power(mrho,4)*(-s + t) + 2*Power(mrho,2)*t*(3*s + t) - 2*Power(mpion,4)*(-5*Power(mrho,2) + s + 12*t) - 
          Power(mpion,2)*(Power(mrho,4) + 5*Power(s,2) - 2*Power(mrho,2)*(s - 6*t) - 12*s*t - 24*Power(t,2))) + 
       Power(eta1,2)*(-4*Power(mpion,6) - Power(s,3) + 2*Power(s,2)*t + 6*s*Power(t,2) + 4*Power(t,3) + 
          Power(mrho,4)*(s + 2*t) - 2*Power(mrho,2)*t*(4*s + 3*t) + 2*Power(mpion,4)*(-5*Power(mrho,2) + s + 6*t) + 
          4*Power(mpion,2)*(Power(s,2) + 4*Power(mrho,2)*t - 2*s*t - 3*Power(t,2))) + 
       Power(eta2,2)*(-4*Power(mpion,6) + 12*Power(mpion,4)*t + 
          Power(mpion,2)*(Power(mrho,4) + Power(s,2) - 4*s*t - 12*Power(t,2) - 2*Power(mrho,2)*(s + 2*t)) + 
          t*(-3*Power(mrho,4) + 2*Power(mrho,2)*(s + 2*t) + Power(s + 2*t,2)))))/
   (16.*Power(mrho,2)*(-Power(mrho,2) + s)*(-Power(ma1,2) + t));

   double MaMd = -(Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(eta2*
         (-4*delta*Power(mpion,6) + 4*Power(mpion,4)*(Power(mrho,2) - 2*C4*Power(mrho,4) + 3*delta*t) + 
           Power(mpion,2)*(-((2 + delta)*Power(mrho,2)*(s + 4*t)) + 2*Power(mrho,4)*(-1 + delta + 8*C4*t) + 
              delta*(Power(s,2) - 4*s*t - 12*Power(t,2))) + 
           t*(-8*C4*Power(mrho,6) + Power(mrho,4)*(6 - 4*delta + 8*C4*(2*s - t)) - 
              Power(mrho,2)*((2 + delta)*s + 8*C4*Power(s,2) - 4*(1 + delta)*t) + 
              delta*(3*Power(s,2) + 4*s*t + 4*Power(t,2)))) + 
        eta1*(4*delta*Power(mpion,6) - 8*C4*Power(mrho,6)*t + 
           2*Power(mrho,4)*(s - delta*s - (-2 + delta)*t + 4*C4*Power(t,2)) + 
           Power(mrho,2)*((-2 + delta)*Power(s,2) + 8*delta*s*t + 8*C4*Power(s,2)*t + 
              2*(-2 + 3*delta)*Power(t,2)) + delta*(Power(s,3) - 4*Power(s,2)*t - 6*s*Power(t,2) - 4*Power(t,3)) - 
           2*Power(mpion,4)*((2 - 5*delta)*Power(mrho,2) - 4*C4*Power(mrho,4) + delta*(s + 6*t)) - 
           2*Power(mpion,2)*(-8*C4*Power(mrho,6) - Power(mrho,2)*((2 + delta)*s + 4*t - 8*delta*t) + 
              2*delta*(Power(s,2) - 2*s*t - 3*Power(t,2)) + Power(mrho,4)*(4 + 8*C4*(s + t))))))/
   (16.*Power(mrho,2)*(-Power(ma1,2) + t));

   double MbMb = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + t,2) - 2*Power(mpion,2)*(Power(mrho,2) + t)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + t,2));

   double MbMc = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(-2*Power(mrho,2) + delta*s)*
     (-4*Power(mpion,4)*(-5*Power(mrho,2) + s) + s*(-3*Power(mrho,4) - s*t + Power(mrho,2)*(s + 3*t + 2*u)) + 
       Power(mpion,2)*(12*Power(mrho,4) + s*(s + 4*t) - Power(mrho,2)*(9*s + 8*t + 12*u))))/
   (16.*Power(mrho,4)*(-Power(mrho,2) + s)*(-Power(mpion,2) + t));

   double MbMd = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(mpion,4)*
        ((-2 + 4*delta)*Power(mrho,2) + 8*C4*Power(mrho,4) + 5*delta*s) - 8*C4*Power(mrho,6)*t + 
       delta*s*t*(s + t) + Power(mrho,2)*(delta*Power(s,2) - (2 + 3*delta)*s*t - 2*Power(t,2)) + 
       2*Power(mrho,4)*((-1 + delta)*s + t + 4*C4*t*(2*s + t)) - 
       Power(mpion,2)*(8*C4*Power(mrho,6) + delta*s*(s + 6*t) + 2*Power(mrho,4)*(-3 + 8*C4*t) + 
          Power(mrho,2)*((-2 + 9*delta)*s + 4*(-1 + delta)*t))))/(16.*Power(mrho,4)*(-Power(mpion,2) + t));
   
   double McMc = (Power(C,2)*Power(ghat,4)*Power(-2*Power(mrho,2) + delta*s,2)*
     (-3*Power(mrho,4)*s + Power(s,3) - 4*Power(mpion,2)*(-3*Power(mrho,4) - 2*Power(mrho,2)*s + Power(s,2)) - 
       2*Power(mrho,2)*(Power(s,2) - Power(t - u,2))))/(16.*Power(mrho,6)*Power(-Power(mrho,2) + s,2));

   double McMd = -(Power(C,2)*Power(ghat,4)*(-2*Power(mrho,2) + delta*s)*
      (-24*C4*Power(mrho,6)*(t - u) + Power(mrho,4)*
         (-6*delta*s + 10*t - delta*t + 8*C4*s*(t - u) - 10*u + delta*u) + 
        delta*s*(2*Power(s,2) - Power(t,2) + Power(u,2)) - 
        2*delta*Power(mpion,2)*(-12*Power(mrho,4) + Power(mrho,2)*(-8*s + t - u) + s*(4*s - t + u)) + 
        Power(mrho,2)*(-4*delta*Power(s,2) + (-2 + delta)*s*(t - u) + delta*(5*Power(t,2) - 8*t*u + 3*Power(u,2)))))
    /(32.*Power(mrho,6)*(-Power(mrho,2) + s));

   double MdMd = (Power(C,2)*Power(ghat,4)*(32*Power(C4,2)*Power(mrho,10) + 4*Power(delta,2)*Power(mpion,4)*s - 
       8*C4*Power(mrho,8)*(6 - 3*delta + 8*C4*s) + 
       Power(mrho,6)*(12 - 8*delta + Power(delta,2) + 32*Power(C4,2)*Power(s,2) - 
          8*C4*(3*(-2 + delta)*s + delta*(5*t - u))) + 2*Power(delta,2)*s*u*(t + u) - 
       2*delta*Power(mrho,4)*((2 + delta)*s - 8*t + delta*t + 2*u - 4*C4*s*(3*t + u)) - 
       2*delta*Power(mpion,2)*(-16*C4*Power(mrho,6) + Power(mrho,4)*(6 - 7*delta + 16*C4*s) + 
          Power(mrho,2)*((2 - 7*delta)*s + delta*(t - u)) + delta*s*(2*s + t + 3*u)) + 
       delta*Power(mrho,2)*(-2*s*(delta*t + 2*(-1 + delta)*u) + delta*(3*Power(t,2) - 4*t*u + Power(u,2)))))/
   (16.*Power(mrho,6));

   double result = MaMa + MaMb + MaMc + MaMd + MbMb + 2*(MbMc + MbMd) 
                   + McMc + 2*McMd + MdMd;
   return(result);
}

/*************************************************************************************/
//rho -> pi + pi + gamma  (C3.1 + C3.2)
/*************************************************************************************/

double C3p1(double s, double t)
{
   double C = 0.059;
   double ghat = g_tilde;
  
   double MaMa = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) + Power(mpion,4)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 4*Power(s,2) - 2*s*t) + 
          Power(s,2)*(-Power(mrho,4) + Power(s,2) - 2*Power(mrho,2)*t + 2*s*t + 2*Power(t,2)) + 
          2*Power(mpion,2)*s*(Power(mrho,4) + Power(mrho,2)*(s + t) - 2*s*(s + t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(Power(mrho,4) + 4*Power(mrho,2)*s + 4*Power(s,2) - 2*s*t) + 
          Power(s,2)*(Power(mrho,4) + Power(s,2) + 2*Power(mrho,2)*(s - t) + 2*s*t + 2*Power(t,2)) - 
          2*Power(mpion,2)*s*(Power(mrho,4) + Power(mrho,2)*(2*s - t) + 2*s*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) - 
          2*Power(mpion,2)*(-Power(mrho,2) + s)*(2*s*(s + t) - Power(mrho,2)*(2*s + t)) + 
          Power(mpion,4)*(3*Power(mrho,4) + 2*s*(2*s - t) + Power(mrho,2)*(-6*s + 2*t)) + 
          s*(-Power(mrho,2) + s)*(Power(s,2) + 2*s*t + 2*Power(t,2) - Power(mrho,2)*(s + 2*t)))))/
   (32.*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

   double MaMb = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (Power(eta2,2)*(Power(mpion,8) - 4*Power(mpion,2)*s*t*(Power(mrho,2) + s + t) + 
          Power(mpion,4)*(2*s*t + Power(mrho,2)*(s + t)) + 
          s*t*(Power(s,2) + 3*s*t + Power(t,2) + Power(mrho,2)*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) + Power(mpion,4)*(2*Power(mrho,4) + 2*s*t - 3*Power(mrho,2)*(s + t)) + 
          s*t*(2*Power(mrho,4) + Power(s,2) + 3*s*t + Power(t,2) - 3*Power(mrho,2)*(s + t)) - 
          2*Power(mpion,2)*(Power(mrho,4)*(s + t) + 2*s*t*(s + t) - 
             Power(mrho,2)*(Power(s,2) + 4*s*t + Power(t,2)))) - 
       2*eta1*eta2*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 2*Power(mpion,4)*s*t + 
          s*t*(Power(s,2) + 3*s*t + Power(t,2) - 2*Power(mrho,2)*(s + t)) + 
          Power(mpion,2)*(-4*s*t*(s + t) + Power(mrho,2)*(Power(s,2) + 4*s*t + Power(t,2))))))/
   (16.*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2))*(-Power(ma1,2) + t));

   double MaMc = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (-(eta2*(Power(mpion,2) + s)*(-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*s) + 
            s*(-Power(mrho,2) + s + 2*t))) + eta1*
        (-4*Power(mpion,6) + s*(-Power(mrho,2) + s)*(-Power(mrho,2) + s + t) + 
          Power(mpion,4)*(3*Power(mrho,2) + s + t) - 
          Power(mpion,2)*(Power(mrho,4) + 2*s*(s - t) + Power(mrho,2)*(-s + t)))))/
   (8.*(-Power(mpion,2) + s)*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

   double MaMd = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (-(eta2*(Power(mpion,2) + s)) + eta1*(-Power(mrho,2) + s + t))*
     (-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*t) + t*(-Power(mrho,2) + 2*s + t)))/
   (8.*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2))*(-Power(mpion,2) + t));

   double MaMe = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(-Power(ma1,2) + s)*
     (eta2*(-2*Power(mpion,4)*(delta - 4*C4*Power(mrho,2))*(Power(mrho,2) + 4*s) + 
          Power(mpion,2)*(-2*Power(mrho,4)*(-2 + delta + 8*C4*s) + 8*delta*s*(s + t) - 
             Power(mrho,2)*((-10 + delta)*s - (-2 + delta)*t + 32*C4*s*(s + t))) + 
          s*(2*Power(mrho,4)*(-2 + delta + 4*C4*s) - 2*delta*Power(s + t,2) + 
             Power(mrho,2)*((-6 + delta)*s + (-2 + delta)*t + 8*C4*Power(s + t,2)))) + 
       eta1*(4*Power(mpion,4)*(6*C4*Power(mrho,4) + 2*delta*s + Power(mrho,2)*(1 - 2*delta - 8*C4*s)) + 
          2*delta*s*Power(s + t,2) - Power(mrho,2)*
           ((-6 + 5*delta)*Power(s,2) + 2*(-2 + 3*delta)*s*t + (-2 + delta)*Power(t,2) + 8*C4*s*Power(s + t,2)) + 
          Power(mrho,4)*((-2 + delta)*(3*s + t) + 8*C4*s*(s + 2*t)) - 
          2*Power(mpion,2)*(4*delta*s*(s + t) - 
             Power(mrho,2)*(-6*s + 7*delta*s - 2*t + 3*delta*t + 16*C4*s*(s + t)) + 
             2*Power(mrho,4)*(-2 + delta + 4*C4*(2*s + t))))))/
   (8.*Power(mrho,2)*(Power(ma1,4) + Power(ma1,2)*(Power(Gammaa1,2) - 2*s) + Power(s,2)));

   double MbMb = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) - Power(mpion,4)*(Power(mrho,4) + 2*Power(mrho,2)*t + 2*(s - 2*t)*t) + 
          Power(t,2)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 2*Power(s,2) + 2*s*t + Power(t,2)) + 
          2*Power(mpion,2)*t*(Power(mrho,4) + Power(mrho,2)*(s + t) - 2*t*(s + t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(Power(mrho,4) + 4*Power(mrho,2)*t - 2*(s - 2*t)*t) + 
          Power(t,2)*(Power(mrho,4) + 2*Power(s,2) - 2*Power(mrho,2)*(s - t) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(Power(mrho,4) - Power(mrho,2)*(s - 2*t) + 2*t*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(3*Power(mrho,4) + 2*Power(mrho,2)*(s - 3*t) - 2*(s - 2*t)*t) + 
          t*(-Power(mrho,2) + t)*(2*Power(s,2) + 2*s*t + Power(t,2) - Power(mrho,2)*(2*s + t)) - 
          2*Power(mpion,2)*(-Power(mrho,2) + t)*(2*t*(s + t) - Power(mrho,2)*(s + 2*t)))))/
   (32.*Power(-Power(ma1,2) + t,2));

   double MbMc = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
     (-(eta2*(Power(mpion,2) + t)) + eta1*(-Power(mrho,2) + s + t))*
     (-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*s) + s*(-Power(mrho,2) + s + 2*t)))/
   (16.*(-Power(mpion,2) + s)*(-Power(ma1,2) + t));

   double MbMd = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
     (-(eta2*(Power(mpion,2) + t)*(-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*t) + 
            t*(-Power(mrho,2) + 2*s + t))) + eta1*
        (-4*Power(mpion,6) + t*(-Power(mrho,2) + t)*(-Power(mrho,2) + s + t) + 
          Power(mpion,4)*(3*Power(mrho,2) + s + t) + 
          Power(mpion,2)*(-Power(mrho,4) + 2*(s - t)*t + Power(mrho,2)*(-s + t)))))/
   (16.*(-Power(ma1,2) + t)*(-Power(mpion,2) + t));

   double MbMe = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(eta2*(-2*Power(mpion,4)*(delta - 4*C4*Power(mrho,2))*
           (Power(mrho,2) + 4*t) + Power(mpion,2)*
           (8*delta*t*(s + t) - 2*Power(mrho,4)*(-2 + delta + 8*C4*t) - 
             Power(mrho,2)*(-((-2 + delta)*s) + (-10 + delta)*t + 32*C4*t*(s + t))) + 
          t*(-2*delta*Power(s + t,2) + 2*Power(mrho,4)*(-2 + delta + 4*C4*t) + 
             Power(mrho,2)*((-2 + delta)*s + (-6 + delta)*t + 8*C4*Power(s + t,2)))) + 
       eta1*(2*delta*t*Power(s + t,2) - Power(mrho,2)*
           ((-2 + delta)*Power(s,2) + 2*(-2 + 3*delta)*s*t + (-6 + 5*delta)*Power(t,2) + 8*C4*t*Power(s + t,2)) + 
          Power(mrho,4)*(8*C4*t*(2*s + t) + (-2 + delta)*(s + 3*t)) + 
          4*Power(mpion,4)*(6*C4*Power(mrho,4) + 2*delta*t + Power(mrho,2)*(1 - 2*delta - 8*C4*t)) - 
          2*Power(mpion,2)*(4*delta*t*(s + t) - 
             Power(mrho,2)*(-2*s + 3*delta*s - 6*t + 7*delta*t + 16*C4*t*(s + t)) + 
             2*Power(mrho,4)*(-2 + delta + 4*C4*(s + 2*t))))))/(16.*Power(mrho,2)*(-Power(ma1,2) + t));

   double McMc = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + s,2) - 2*Power(mpion,2)*(Power(mrho,2) + s)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + s,2));

   double McMd = (Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*(Power(mrho,6) - 2*Power(mrho,4)*(s + t) + s*t*(s + t) + 
       Power(mpion,4)*(-Power(mrho,2) + s + t) + Power(mrho,2)*(Power(s,2) + s*t + Power(t,2)) - 
       Power(mpion,2)*(2*Power(mrho,4) - 3*Power(mrho,2)*(s + t) + Power(s + t,2))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + s)*(-Power(mpion,2) + t));

   double McMe = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(mpion,4)*(2 - 3*delta + 8*C4*Power(mrho,2)) + 
       Power(mrho,4)*(-2 + delta + 8*C4*s) + Power(mpion,2)*
        (-8*C4*Power(mrho,4) - (2 + 3*delta)*s + 4*Power(mrho,2)*(1 + 4*C4*s) + (-2 + delta)*t) + 
       s*(2*delta*s + (2 + 3*delta)*t) - Power(mrho,2)*((-2 + 3*delta)*s + (-2 + delta)*t + 8*C4*s*(s + 2*t))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + s));

   double MdMd = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + t,2) - 2*Power(mpion,2)*(Power(mrho,2) + t)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + t,2));

   double MdMe = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(Power(mpion,4)*(2 - 3*delta + 8*C4*Power(mrho,2)) + 
       Power(mrho,4)*(-2 + delta + 8*C4*t) + t*((2 + 3*delta)*s + 2*delta*t) - 
       Power(mrho,2)*((-2 + delta)*s + (-2 + 3*delta)*t + 8*C4*t*(2*s + t)) + 
       Power(mpion,2)*(-8*C4*Power(mrho,4) + (-2 + delta)*s - (2 + 3*delta)*t + 4*Power(mrho,2)*(1 + 4*C4*t))))/
   (8.*Power(mrho,2)*(-Power(mpion,2) + t));

   double MeMe = (Power(C,2)*Power(ghat,4)*(8*Power(mpion,4)*Power(delta - 4*C4*Power(mrho,2),2) + 
       2*Power(delta,2)*Power(s + t,2) - 2*delta*Power(mrho,2)*(s + t)*(3*(-2 + delta) + 8*C4*(s + t)) + 
       Power(mrho,4)*(3*Power(-2 + delta,2) + 24*C4*(-2 + delta)*(s + t) + 32*Power(C4,2)*Power(s + t,2)) - 
       4*Power(mpion,2)*(delta - 4*C4*Power(mrho,2))*
        (2*delta*(s + t) - Power(mrho,2)*(3*(-2 + delta) + 8*C4*(s + t)))))/(4.*Power(mrho,4));

   double result = MaMa + MaMb + MaMc + MaMd + MaMe + MbMb + 2*(MbMc + MbMd + MbMe) 
                   + McMc + 2*(McMd + McMe) + MdMd + 2*MdMe + MeMe;
   return(result);

}

double C3p2(double s, double t)
{
   double C = 0.059;
   double ghat = g_tilde;
  
   double u = -s - t + 2*mpion*mpion + mrho*mrho; 

   double MaMa = (Power(C,2)*Power(eta1 - eta2,2)*Power(ghat,4)*
     (-2*eta1*eta2*(Power(mpion,8) - Power(mpion,4)*(Power(mrho,4) + 2*Power(mrho,2)*t + 2*(s - 2*t)*t) + 
          Power(t,2)*(-Power(mrho,4) - 2*Power(mrho,2)*s + 2*Power(s,2) + 2*s*t + Power(t,2)) + 
          2*Power(mpion,2)*t*(Power(mrho,4) + Power(mrho,2)*(s + t) - 2*t*(s + t))) + 
       Power(eta2,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(Power(mrho,4) + 4*Power(mrho,2)*t - 2*(s - 2*t)*t) + 
          Power(t,2)*(Power(mrho,4) + 2*Power(s,2) - 2*Power(mrho,2)*(s - t) + 2*s*t + Power(t,2)) - 
          2*Power(mpion,2)*t*(Power(mrho,4) - Power(mrho,2)*(s - 2*t) + 2*t*(s + t))) + 
       Power(eta1,2)*(Power(mpion,8) - 2*Power(mpion,6)*Power(mrho,2) + 
          Power(mpion,4)*(3*Power(mrho,4) + 2*Power(mrho,2)*(s - 3*t) - 2*(s - 2*t)*t) + 
          t*(-Power(mrho,2) + t)*(2*Power(s,2) + 2*s*t + Power(t,2) - Power(mrho,2)*(2*s + t)) - 
          2*Power(mpion,2)*(-Power(mrho,2) + t)*(2*t*(s + t) - Power(mrho,2)*(s + 2*t)))))/
   (32.*Power(-Power(ma1,2) + t,2));

   double MaMb = (Power(C,2)*(-2 + delta)*(eta1 - eta2)*Power(ghat,4)*
     (-(eta2*(Power(mpion,2) + t)*(-Power(mpion,4) + Power(mpion,2)*(Power(mrho,2) - 2*t) + 
            t*(-Power(mrho,2) + 2*s + t))) + eta1*
        (-4*Power(mpion,6) + t*(-Power(mrho,2) + t)*(-Power(mrho,2) + s + t) + 
          Power(mpion,4)*(3*Power(mrho,2) + s + t) + 
          Power(mpion,2)*(-Power(mrho,4) + 2*(s - t)*t + Power(mrho,2)*(-s + t)))))/
   (8.*(-Power(ma1,2) + t)*(-Power(mpion,2) + t));

   double MaMc = -(Power(C,2)*Power(ghat,4)*(-2*delta*Power(mpion,2) - (-2 + delta)*Power(mrho,2) + delta*(s + t))*
      (Power(eta2,2)*(-4*Power(mpion,4)*(s - t) + (s - t)*(-4*Power(mrho,2) + s - t)*t + 
           Power(mpion,2)*(Power(s,2) + 2*s*t - 3*Power(t,2))) + 
        Power(eta1,2)*(8*Power(mpion,6) + Power(s,3) + 2*Power(mrho,4)*(s - t) + 5*Power(s,2)*t + s*Power(t,2) + 
           Power(t,3) - 2*Power(mpion,2)*(s + t)*(-2*Power(mrho,2) + s + t) - 
           2*Power(mpion,4)*(2*Power(mrho,2) + 3*s + t) + Power(mrho,2)*(-3*Power(s,2) - 2*s*t + Power(t,2))) + 
        eta1*eta2*(-8*Power(mpion,6) - Power(s,3) - 2*Power(mrho,4)*(s - t) + 
           2*Power(mpion,4)*(2*Power(mrho,2) + 5*s - t) - 6*Power(s,2)*t + s*Power(t,2) - 2*Power(t,3) + 
           Power(mrho,2)*(3*Power(s,2) + 6*s*t - 5*Power(t,2)) + 
           Power(mpion,2)*(Power(s,2) + 2*s*t + 5*Power(t,2) - 4*Power(mrho,2)*(s + t)))))/
   (16.*Power(mrho,2)*(-Power(ma1,2) + t)*(-2*Power(mpion,2) + s + t));

   double MaMd = (Power(C,2)*(eta1 - eta2)*Power(ghat,4)*(-(eta2*
          (-2*Power(mpion,4)*(4*C4*Power(mrho,4) + 2*delta*(s - 3*t) - Power(mrho,2)*(delta - 16*C4*t)) + 
            Power(mpion,2)*(2*Power(mrho,4)*(-2 + delta + 8*C4*t) + delta*(Power(s,2) - 6*s*t - 11*Power(t,2)) + 
               Power(mrho,2)*(-((-2 + delta)*s) + (-10 + delta)*t + 32*C4*t*(s + t))) + 
            t*(-2*Power(mrho,4)*(-2 + delta + 4*C4*t) + delta*(3*Power(s,2) + 2*s*t + 3*Power(t,2)) + 
               Power(mrho,2)*((2 - 5*delta)*s + 3*(2 + delta)*t - 8*C4*Power(s + t,2))))) + 
       eta1*(8*delta*Power(mpion,6) + delta*(Power(s,3) + 7*Power(s,2)*t + 5*s*Power(t,2) + 3*Power(t,3)) - 
          2*Power(mrho,2)*((-1 + 2*delta)*Power(s,2) + 2*(-1 + 2*delta)*s*t + (-3 + 2*delta)*Power(t,2) + 
             4*C4*t*Power(s + t,2)) + Power(mrho,4)*((-2 + 3*delta)*s + (-6 + delta)*t + 8*C4*t*(2*s + t)) + 
          Power(mpion,4)*(24*C4*Power(mrho,4) + 6*delta*(-s + t) - 4*Power(mrho,2)*(-1 + 3*delta + 8*C4*t)) - 
          2*Power(mpion,2)*(delta*(Power(s,2) + 6*s*t + 5*Power(t,2)) - 
             Power(mrho,2)*(-2*s + 5*delta*s - 6*t + 9*delta*t + 16*C4*t*(s + t)) + 
             2*Power(mrho,4)*(-2 + delta + 4*C4*(s + 2*t))))))/(16.*Power(mrho,2)*(-Power(ma1,2) + t));

   double MbMb = -(Power(C,2)*Power(-2 + delta,2)*Power(ghat,4)*Power(mpion,2)*
      (Power(mpion,4) + Power(-Power(mrho,2) + t,2) - 2*Power(mpion,2)*(Power(mrho,2) + t)))/
   (4.*Power(mrho,2)*Power(-Power(mpion,2) + t,2));

   double MbMc = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(-2*Power(mrho,2) + delta*u)*
     (-4*Power(mpion,4)*(-5*Power(mrho,2) + u) + u*(-3*Power(mrho,4) - t*u + Power(mrho,2)*(2*s + 3*t + u)) + 
       Power(mpion,2)*(12*Power(mrho,4) + u*(4*t + u) - Power(mrho,2)*(12*s + 8*t + 9*u))))/
   (16.*Power(mrho,4)*(-Power(mpion,2) + t)*(-Power(mrho,2) + u));

   double MbMd = (Power(C,2)*(-2 + delta)*Power(ghat,4)*(6*delta*Power(mpion,6) + delta*s*t*(s + t) + 
       Power(mrho,6)*(-2 + 3*delta + 8*C4*t) + 
       Power(mrho,2)*(delta*Power(s,2) + (2 + 3*delta)*s*t + 3*delta*Power(t,2)) - 
       2*Power(mrho,4)*(-s + 2*delta*s - t + 3*delta*t + 4*C4*t*(2*s + t)) - 
       Power(mpion,4)*((-2 + 9*delta)*Power(mrho,2) - 8*C4*Power(mrho,4) + delta*(s + 9*t)) - 
       Power(mpion,2)*(8*C4*Power(mrho,6) + 2*Power(mrho,4)*(-2 + delta - 8*C4*t) + 
          Power(mrho,2)*(2*s - 7*delta*s + 2*t + 5*delta*t) + delta*(Power(s,2) - 3*Power(t,2)))))/
   (16.*Power(mrho,4)*(-Power(mpion,2) + t));

   double McMc = (Power(C,2)*Power(ghat,4)*Power(-2*Power(mrho,2) + delta*u,2)*
     (-3*Power(mrho,4)*u + Power(u,3) + 2*Power(mrho,2)*(Power(s,2) - 2*s*t + Power(t,2) - Power(u,2)) - 
       4*Power(mpion,2)*(-3*Power(mrho,4) - 2*Power(mrho,2)*u + Power(u,2))))/
   (16.*Power(mrho,6)*Power(-Power(mrho,2) + u,2));

   double McMd = -(Power(C,2)*Power(ghat,4)*(-2*Power(mrho,2) + delta*u)*
      (24*C4*Power(mrho,6)*(s - t) + Power(mrho,4)*
         ((-10 + delta)*s + 10*t - delta*t - 6*delta*u + 8*C4*(-s + t)*u) + 
        delta*u*(Power(s,2) - Power(t,2) + 2*Power(u,2)) - 
        2*delta*Power(mpion,2)*(-12*Power(mrho,4) + Power(mrho,2)*(-s + t - 8*u) + u*(s - t + 4*u)) + 
        Power(mrho,2)*(3*delta*Power(s,2) + 5*delta*Power(t,2) + (-2 + delta)*t*u - 4*delta*Power(u,2) - 
           s*(8*delta*t + (-2 + delta)*u))))/(32.*Power(mrho,6)*(-Power(mrho,2) + u));

   double MdMd = (Power(C,2)*Power(ghat,4)*(32*Power(C4,2)*Power(mrho,10) + 4*Power(delta,2)*Power(mpion,4)*u + 
       2*Power(delta,2)*s*(s + t)*u - 8*C4*Power(mrho,8)*(6 - 3*delta + 8*C4*u) - 
       2*delta*Power(mrho,4)*(2*s - 8*t + delta*t + 2*u + delta*u - 4*C4*(s + 3*t)*u) + 
       delta*Power(mrho,2)*(delta*Power(s,2) + delta*t*(3*t - 2*u) - 4*s*(delta*t + (-1 + delta)*u)) - 
       2*delta*Power(mpion,2)*(-16*C4*Power(mrho,6) + delta*u*(3*s + t + 2*u) + 
          Power(mrho,4)*(6 - 7*delta + 16*C4*u) + Power(mrho,2)*(-(delta*s) + delta*t + 2*u - 7*delta*u)) + 
       Power(mrho,6)*(12 - 8*delta + Power(delta,2) + 32*Power(C4,2)*Power(u,2) - 
          8*C4*(-(delta*s) + 5*delta*t - 6*u + 3*delta*u))))/(16.*Power(mrho,6));

   double result = MaMa + MaMb + MaMc + MaMd + MbMb + 2*(MbMc + MbMd) 
                   + McMc + 2*(McMd) + MdMd;
   return(result);

}


/*************************************************************************************/
//pi + Kstar -> K + gamma  (C4.1 + C4.2)
/*************************************************************************************/

double C4p1(double s, double t)
{
   double C = 0.059;
   double gk = 9.2;

   double MaMa = -(Power(C,2)*Power(gk,4)*Power(mpion,2)*(Power(mK,4) + Power(-Power(mKstar,2) + t,2) - 2*Power(mK,2)*(Power(mKstar,2) + t)))/
   (2.*Power(mKstar,2)*Power(-Power(mpion,2) + t,2));

   double MaMb = (Power(C,2)*Power(gk,4)*(2*Power(mKstar,6) - 3*Power(mKstar,4)*s + Power(mKstar,2)*Power(s,2) - 3*Power(mKstar,4)*t + Power(mKstar,2)*s*t + Power(s,2)*t + 
       2*s*Power(t,2) + Power(t,3) - Power(mpion,2)*(Power(mKstar,2) + s - t)*(-Power(mKstar,2) + t) + 
       Power(mK,4)*(Power(mKstar,2) + 4*Power(mpion,2) - s + t) - 
       Power(mK,2)*(3*Power(mKstar,4) + Power(s,2) + s*t + 2*Power(t,2) - Power(mKstar,2)*(4*s + t) + Power(mpion,2)*(7*Power(mKstar,2) - s + 5*t))))/
   (16.*Power(mKstar,2)*(-Power(mpion,2) + t)*(-Power(mK,2) - Power(mpion,2) + s + t));

   double MaMc = (Power(C,2)*Power(gk,4)*(Power(mKstar,6) - 2*Power(mKstar,4)*s + Power(mKstar,2)*Power(s,2) - 2*Power(mKstar,4)*t + Power(mKstar,2)*s*t + Power(s,2)*t + 
       Power(mKstar,2)*Power(t,2) + s*Power(t,2) - Power(mK,2)*(-Power(mKstar,2) - Power(mpion,2) + s)*(-Power(mKstar,2) + s + t) - 
       Power(mpion,2)*(-Power(mKstar,2) + t)*(-Power(mKstar,2) + s + t)))/(8.*Power(mKstar,2)*(-Power(mK,2) + s)*(-Power(mpion,2) + t));

   double MaMd = (-3*Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + Power(mK,2)*(Power(mKstar,2) + Power(mpion,2) - s) + Power(mKstar,2)*s + 
       Power(mpion,2)*(Power(mKstar,2) - t) + Power(mKstar,2)*t + s*t))/(16.*Power(mKstar,2)*(-Power(mpion,2) + t));

   double MbMb = -(Power(C,2)*Power(gk,4)*(-Power(mK,6) + 4*Power(mKstar,6) - Power(mpion,6) - 4*Power(mKstar,4)*s - 3*Power(mKstar,2)*Power(s,2) + Power(s,3) + 
        Power(mpion,4)*(Power(mKstar,2) + 3*s - t) - 4*Power(mKstar,4)*t + 2*Power(mKstar,2)*s*t + 3*Power(s,2)*t - 3*Power(mKstar,2)*Power(t,2) + 
        3*s*Power(t,2) + Power(t,3) + Power(mK,4)*(Power(mKstar,2) + 5*Power(mpion,2) - s + 3*t) + 
        Power(mpion,2)*(-4*Power(mKstar,4) - 3*Power(s,2) - 2*Power(mKstar,2)*(s - 3*t) - 2*s*t + Power(t,2)) + 
        Power(mK,2)*(-4*Power(mKstar,4) + 5*Power(mpion,4) + Power(s,2) + Power(mKstar,2)*(6*s - 2*t) - 2*s*t - 3*Power(t,2) - 
           6*Power(mpion,2)*(Power(mKstar,2) + s + t))))/(32.*Power(mKstar,2)*Power(-Power(mK,2) - Power(mpion,2) + s + t,2));

   double MbMc = -(Power(C,2)*Power(gk,4)*(2*Power(mKstar,6) - 3*Power(mKstar,4)*s + Power(s,3) + Power(mpion,4)*(Power(mKstar,2) + s - t) - 3*Power(mKstar,4)*t + 
        Power(mKstar,2)*s*t + 2*Power(s,2)*t + Power(mKstar,2)*Power(t,2) + s*Power(t,2) + 
        Power(mK,2)*(4*Power(mpion,4) + (-Power(mKstar,2) + s)*(-Power(mKstar,2) + s - t) + Power(mpion,2)*(-7*Power(mKstar,2) - 5*s + t)) - 
        Power(mpion,2)*(3*Power(mKstar,4) + 2*Power(s,2) + s*t + Power(t,2) - Power(mKstar,2)*(s + 4*t))))/
   (32.*Power(mKstar,2)*(-Power(mK,2) + s)*(-Power(mK,2) - Power(mpion,2) + s + t));

   double MbMd = (3*Power(C,2)*Power(gk,4)*(-Power(mK,4) + Power(mpion,4) + 4*Power(mKstar,2)*s + Power(s,2) - 2*Power(mpion,2)*(-2*Power(mKstar,2) + s) - 
       4*Power(mKstar,2)*t - Power(t,2) + 2*Power(mK,2)*(-2*Power(mKstar,2) + t)))/(64.*Power(mKstar,2)*(-Power(mK,2) - Power(mpion,2) + s + t));

   double McMc = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mpion,4) + Power(-Power(mKstar,2) + s,2) - 2*Power(mpion,2)*(Power(mKstar,2) + s)))/
   (8.*Power(mKstar,2)*Power(-Power(mK,2) + s,2));

   double McMd = (-3*Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + Power(mK,2)*(Power(mKstar,2) + Power(mpion,2) - s) + Power(mKstar,2)*s + 
       Power(mpion,2)*(Power(mKstar,2) - t) + Power(mKstar,2)*t + s*t))/(32.*Power(mKstar,2)*(-Power(mK,2) + s));

   double MdMd = (27*Power(C,2)*Power(gk,4))/32.;

   double result = MaMa + 2*(MaMb + MaMc + MaMd) + MbMb + 2*(MbMc + MbMd) + McMc 
                   + 2*McMd + MdMd;

   return(result);
}

double C4p2(double s, double t)
{
   double C = 0.059;
   double gk = 9.2;

   double MaMa = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mpion,4) + Power(-Power(mKstar,2) + s,2) - 2*Power(mpion,2)*(Power(mKstar,2) + s)))/
   (16.*Power(mKstar,2)*Power(-Power(mK,2) + s,2));

   double MaMb = (Power(C,2)*Power(gk,4)*(2*Power(mKstar,6) - 3*Power(mKstar,4)*s + Power(s,3) + Power(mpion,4)*(Power(mKstar,2) + s - t) - 3*Power(mKstar,4)*t + 
       Power(mKstar,2)*s*t + 2*Power(s,2)*t + Power(mKstar,2)*Power(t,2) + s*Power(t,2) + 
       Power(mK,2)*(4*Power(mpion,4) + (-Power(mKstar,2) + s)*(-Power(mKstar,2) + s - t) + Power(mpion,2)*(-7*Power(mKstar,2) - 5*s + t)) - 
       Power(mpion,2)*(3*Power(mKstar,4) + 2*Power(s,2) + s*t + Power(t,2) - Power(mKstar,2)*(s + 4*t))))/
   (64.*Power(mKstar,2)*(-Power(mK,2) + s)*(-Power(mK,2) - Power(mpion,2) + s + t));

   double MaMc = -(Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + Power(mK,2)*(Power(mKstar,2) + Power(mpion,2) - s) + Power(mKstar,2)*s + Power(mpion,2)*(Power(mKstar,2) - t) + 
        Power(mKstar,2)*t + s*t))/(64.*Power(mKstar,2)*(-Power(mK,2) + s));

   double MbMb = -(Power(C,2)*Power(gk,4)*(-Power(mK,6) + 4*Power(mKstar,6) - Power(mpion,6) - 4*Power(mKstar,4)*s - 3*Power(mKstar,2)*Power(s,2) + Power(s,3) + 
        Power(mpion,4)*(Power(mKstar,2) + 3*s - t) - 4*Power(mKstar,4)*t + 2*Power(mKstar,2)*s*t + 3*Power(s,2)*t - 3*Power(mKstar,2)*Power(t,2) + 
        3*s*Power(t,2) + Power(t,3) + Power(mK,4)*(Power(mKstar,2) + 5*Power(mpion,2) - s + 3*t) + 
        Power(mpion,2)*(-4*Power(mKstar,4) - 3*Power(s,2) - 2*Power(mKstar,2)*(s - 3*t) - 2*s*t + Power(t,2)) + 
        Power(mK,2)*(-4*Power(mKstar,4) + 5*Power(mpion,4) + Power(s,2) + Power(mKstar,2)*(6*s - 2*t) - 2*s*t - 3*Power(t,2) - 
           6*Power(mpion,2)*(Power(mKstar,2) + s + t))))/(64.*Power(mKstar,2)*Power(-Power(mK,2) - Power(mpion,2) + s + t,2));

   double MbMc = (Power(C,2)*Power(gk,4)*(Power(mK,4) - Power(mpion,4) - 4*Power(mKstar,2)*s - Power(s,2) + 2*Power(mpion,2)*(-2*Power(mKstar,2) + s) + 4*Power(mKstar,2)*t + 
       Power(t,2) - 2*Power(mK,2)*(-2*Power(mKstar,2) + t)))/(128.*Power(mKstar,2)*(-Power(mK,2) - Power(mpion,2) + s + t));

   double McMc = (3*Power(C,2)*Power(gk,4))/64.;

   double result = MaMa + 2*(MaMb + MaMc) + MbMb + 2*MbMc + McMc;

   return(result);
}


/*************************************************************************************/
//pi + K -> Kstar + gamma  (C5.1 + C5.2)
/*************************************************************************************/

double C5p1(double s, double t)
{
  double C = 0.059;
  double gk = 9.2;
  double u = -s - t + mK*mK + mpion*mpion + mKstar*mKstar;

  double MaMa = -(Power(C,2)*Power(gk,4)*Power(mpion,2)*(Power(mK,4) + Power(-Power(mKstar,2) + t,2) - 
        2*Power(mK,2)*(Power(mKstar,2) + t)))/(2.*Power(mKstar,2)*Power(-Power(mpion,2) + t,2));

  double MaMb = -(Power(C,2)*Power(gk,4)*(-2*Power(mK,6) + 2*Power(mKstar,2)*Power(mpion,4) + 
        Power(mK,4)*(3*Power(mKstar,2) + 2*Power(mpion,2) + 3*s + 4*t) + 
        s*(Power(mKstar,4) + Power(mKstar,2)*(s - t) + s*t) + 
        Power(mpion,2)*(Power(mKstar,4) - 3*Power(mKstar,2)*(s + t) + t*(-s + 2*t)) - 
        Power(mK,2)*(Power(mKstar,4) + Power(s,2) + 3*s*t + 2*Power(t,2) + Power(mKstar,2)*(4*s + t) + 
           Power(mpion,2)*(Power(mKstar,2) - s + 4*t))))/
   (16.*Power(mKstar,2)*(-Power(mKstar,2) + s)*(-Power(mpion,2) + t));

   double MaMc = -(Power(C,2)*Power(gk,4)*(Power(mK,2) + Power(mpion,2) - s)*
      (-Power(mKstar,4) + 2*Power(mKstar,2)*s + Power(mKstar,2)*t + Power(mpion,2)*(-3*Power(mKstar,2) + t) + 
        Power(mKstar,2)*u - t*u + Power(mK,2)*(-3*Power(mKstar,2) - Power(mpion,2) + u)))/
   (8.*Power(mKstar,2)*(-Power(mpion,2) + t)*(-Power(mK,2) + u));

   double MaMd = (3*Power(C,2)*Power(gk,4)*(Power(mK,4) + Power(mKstar,2)*(-2*Power(mpion,2) + s - t) + t*(s + t) - 
       Power(mK,2)*(Power(mKstar,2) + s + 2*t)))/(16.*Power(mKstar,2)*(-Power(mpion,2) + t));

   double MbMb = (Power(C,2)*Power(gk,4)*(-3*Power(mKstar,4)*s - 2*Power(mKstar,2)*Power(s,2) + Power(s,3) + 
       2*Power(mK,4)*(-2*Power(mKstar,2) + s) + 2*Power(mpion,4)*(-2*Power(mKstar,2) + s) + 
       2*Power(mKstar,2)*Power(t,2) - 2*Power(mK,2)*
        (2*Power(mpion,2)*(-2*Power(mKstar,2) + s) + (-3*Power(mKstar,2) + s)*(Power(mKstar,2) + s + t - u)) - 
       4*Power(mKstar,2)*t*u + 2*Power(mKstar,2)*Power(u,2) - 
       2*Power(mpion,2)*(-3*Power(mKstar,2) + s)*(Power(mKstar,2) + s - t + u)))/
   (32.*Power(Power(mKstar,3) - mKstar*s,2));

   double MbMc = (Power(C,2)*Power(gk,4)*(Power(mpion,4)*(Power(mKstar,2) + s - 2*t) + 
       2*Power(mK,4)*(-3*Power(mKstar,2) - Power(mpion,2) + u) + 
       Power(mpion,2)*(-5*Power(mKstar,4) - Power(s,2) - s*u + 2*t*u + Power(mKstar,2)*(2*s + 4*t + u)) + 
       s*(3*Power(mKstar,4) + s*u - Power(mKstar,2)*(s + 2*t + 3*u)) + 
       Power(mK,2)*(-7*Power(mKstar,4) + 2*Power(mpion,4) + 
          Power(mpion,2)*(-15*Power(mKstar,2) + 3*s + 2*t - 2*u) - (3*s + 2*t)*u + Power(mKstar,2)*(7*s + 8*t + 7*u)
          )))/(32.*Power(mKstar,2)*(-Power(mKstar,2) + s)*(-Power(mK,2) + u));

   double MbMd = (3*Power(C,2)*Power(gk,4)*(Power(mK,2)*(3*Power(mKstar,2) + s) - Power(mpion,2)*(3*Power(mKstar,2) + s) - 
       (-5*Power(mKstar,2) + s)*(t - u)))/(64.*Power(mKstar,2)*(-Power(mKstar,2) + s));

   double McMc = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mpion,4) + Power(-Power(mKstar,2) + u,2) - 
        2*Power(mpion,2)*(Power(mKstar,2) + u)))/(8.*Power(mKstar,2)*Power(-Power(mK,2) + u,2));

   double McMd = (3*Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + 2*Power(mKstar,2)*s + Power(mKstar,2)*t + 
       Power(mpion,2)*(-3*Power(mKstar,2) + t) + Power(mKstar,2)*u - t*u + 
       Power(mK,2)*(-3*Power(mKstar,2) - Power(mpion,2) + u)))/(32.*Power(mKstar,2)*(-Power(mK,2) + u));

   double MdMd = (27*Power(C,2)*Power(gk,4))/32.;

   double result = MaMa + 2*(MaMb + MaMc + MaMd) + MbMb + 2*(MbMc + MbMd) + McMc + 2*McMd + MdMd;

   return(result);
}

double C5p2(double s, double t)
{
  double C = 0.059;
  double gk = 9.2;
  double u = -s - t + mK*mK + mpion*mpion + mKstar*mKstar;

  double MaMa = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mpion,4) + Power(-Power(mKstar,2) + u,2) - 
        2*Power(mpion,2)*(Power(mKstar,2) + u)))/(16.*Power(mKstar,2)*Power(-Power(mK,2) + u,2));

  double MaMb = -(Power(C,2)*Power(gk,4)*(-2*Power(mK,6) - 2*Power(mKstar,2)*Power(mpion,4) + 
        Power(mK,4)*(-3*Power(mKstar,2) + 2*Power(mpion,2) + 5*s + 4*t) + 
        Power(mpion,2)*(4*Power(mKstar,4) - Power(s,2) + Power(mKstar,2)*(3*s - 5*t) + s*t + 2*Power(t,2)) + 
        s*(s*(s + t) - Power(mKstar,2)*(3*s + t)) - 
        Power(mK,2)*(4*Power(s,2) + 5*s*t + 2*Power(t,2) - Power(mKstar,2)*(6*s + t) + 
           Power(mpion,2)*(-9*Power(mKstar,2) + s + 4*t))))/
   (64.*Power(mKstar,2)*(-Power(mKstar,2) + s)*(-Power(mKstar,2) - Power(mpion,2) + s + t));

   double MaMc = (Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + 2*Power(mKstar,2)*s + Power(mKstar,2)*t + 
       Power(mpion,2)*(-3*Power(mKstar,2) + t) + Power(mKstar,2)*u - t*u + 
       Power(mK,2)*(-3*Power(mKstar,2) - Power(mpion,2) + u)))/(64.*Power(mKstar,2)*(-Power(mK,2) + u));

   double MbMb = (Power(C,2)*Power(gk,4)*(2*Power(mKstar,6) + 4*Power(mKstar,2)*Power(mpion,4) - 7*Power(mKstar,4)*s + Power(s,3) + 
       4*Power(mK,4)*(-2*Power(mKstar,2) + s) - 8*Power(mKstar,4)*t + 8*Power(mKstar,2)*s*t + 
       8*Power(mKstar,2)*Power(t,2) - 4*Power(mK,2)*
        (-Power(mKstar,4) + Power(mpion,2)*(-3*Power(mKstar,2) + s) + s*(s + t) - Power(mKstar,2)*(2*s + t)) + 
       4*Power(mpion,2)*(4*Power(mKstar,4) + s*t - Power(mKstar,2)*(2*s + 5*t))))/
   (64.*Power(Power(mKstar,3) - mKstar*s,2));

   double MbMc = (Power(C,2)*Power(gk,4)*(-(Power(mK,2)*(3*Power(mKstar,2) + s)) + Power(mpion,2)*(3*Power(mKstar,2) + s) + 
       (-5*Power(mKstar,2) + s)*(t - u)))/(128.*Power(mKstar,2)*(-Power(mKstar,2) + s));

   double McMc = (3*Power(C,2)*Power(gk,4))/64.;

   double result = MaMa + 2*(MaMb + MaMc) + MbMb + 2*MbMc + McMc;

   return(result);
}

/*************************************************************************************/
//rho + K -> K + gamma  (C6.1 + C6.2)
/*************************************************************************************/

double C6p1(double s, double t)
{
  double C = 0.059;
  double gk = 9.2;

  double MaMa = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mK,4) + Power(-Power(mrho,2) + s,2) - 
        2*Power(mK,2)*(Power(mrho,2) + s)))/(8.*Power(mrho,2)*Power(-Power(mK,2) + s,2));

  double MaMb = (Power(C,2)*Power(gk,4)*(4*Power(mK,6) + 2*Power(mrho,6) - 2*Power(mK,4)*(3*Power(mrho,2) + 2*s) - 
       3*Power(mrho,4)*(s + t) + Power(mrho,2)*t*(s + t) + s*Power(s + t,2) - 
       Power(mK,2)*(2*Power(mrho,4) + Power(mrho,2)*(s - 5*t) + Power(s + t,2))))/
       (16.*Power(mrho,2)*(-Power(mK,2) + s)*(-2*Power(mK,2) + s + t));

  double MaMc = -(Power(C,2)*Power(gk,4)*(Power(mrho,6) - 2*Power(mrho,4)*(s + t) + s*t*(s + t) + 
        Power(mK,4)*(-Power(mrho,2) + s + t) + Power(mrho,2)*(Power(s,2) + s*t + Power(t,2)) - 
        Power(mK,2)*(2*Power(mrho,4) - 3*Power(mrho,2)*(s + t) + Power(s + t,2))))/
   (16.*Power(mrho,2)*(-Power(mK,2) + s)*(-Power(mK,2) + t));

  double MbMb = -(Power(C,2)*Power(gk,4)*(8*Power(mK,6) + 4*Power(mrho,6) - 4*Power(mrho,4)*(s + t) + Power(s + t,3) - 
        4*Power(mK,4)*(Power(mrho,2) + s + t) + Power(mrho,2)*(-3*Power(s,2) + 2*s*t - 3*Power(t,2)) - 
        2*Power(mK,2)*(4*Power(mrho,4) - 2*Power(mrho,2)*(s + t) + Power(s + t,2))))/
   (8.*Power(mrho,2)*Power(-2*Power(mK,2) + s + t,2));

  double MbMc = (Power(C,2)*Power(gk,4)*(4*Power(mK,6) + 2*Power(mrho,6) - 3*Power(mrho,4)*(s + t) + Power(mrho,2)*s*(s + t) + 
       t*Power(s + t,2) - 2*Power(mK,4)*(3*Power(mrho,2) + 2*t) - 
       Power(mK,2)*(2*Power(mrho,4) + Power(mrho,2)*(-5*s + t) + Power(s + t,2))))/
   (16.*Power(mrho,2)*(-Power(mK,2) + t)*(-2*Power(mK,2) + s + t));

  double McMc = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mK,4) + Power(-Power(mrho,2) + t,2) - 
        2*Power(mK,2)*(Power(mrho,2) + t)))/(8.*Power(mrho,2)*Power(-Power(mK,2) + t,2));

  double result = MaMa + 2*(MaMb + MaMc) + MbMb + 2*MbMc + McMc;

  return(result);
}

double C6p2(double s, double t)
{
  double C = 0.059;
  double gk = 9.2;

  double MaMa = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mK,4) + Power(-Power(mrho,2) + s,2) - 
        2*Power(mK,2)*(Power(mrho,2) + s)))/(16.*Power(mrho,2)*Power(-Power(mK,2) + s,2));

  double MaMb = (Power(C,2)*Power(gk,4)*(Power(mrho,6) - 2*Power(mrho,4)*(s + t) + s*t*(s + t) + 
       Power(mK,4)*(-Power(mrho,2) + s + t) + Power(mrho,2)*(Power(s,2) + s*t + Power(t,2)) - 
       Power(mK,2)*(2*Power(mrho,4) - 3*Power(mrho,2)*(s + t) + Power(s + t,2))))/
   (32.*Power(mrho,2)*(-Power(mK,2) + s)*(-Power(mK,2) + t));

  double MaMc = -(Power(C,2)*Power(gk,4)*(Power(mK,4) - Power(mrho,4) + s*t + Power(mrho,2)*(s + t) - 
        Power(mK,2)*(-2*Power(mrho,2) + s + t)))/(32.*Power(mrho,2)*(-Power(mK,2) + s));

  double MbMb = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mK,4) + Power(-Power(mrho,2) + t,2) - 
        2*Power(mK,2)*(Power(mrho,2) + t)))/(16.*Power(mrho,2)*Power(-Power(mK,2) + t,2));

  double MbMc = -(Power(C,2)*Power(gk,4)*(Power(mK,4) - Power(mrho,4) + s*t + Power(mrho,2)*(s + t) - 
        Power(mK,2)*(-2*Power(mrho,2) + s + t)))/(32.*Power(mrho,2)*(-Power(mK,2) + t));

  double McMc = (3*Power(C,2)*Power(gk,4))/16.;

  double result = MaMa + 2*(MaMb + MaMc) + MbMb + 2*MbMc + McMc;
  return(result);
}


/*************************************************************************************/
//K + Kstar -> pi + gamma  (C7.1 + C7.2)
/*************************************************************************************/

double C7p1(double s, double t)
{
  double C = 0.059;
  double gk = 9.2;

  double MaMa = -(Power(C,2)*Power(gk,4)*Power(mpion,2)*(Power(mK,4) + Power(-Power(mKstar,2) + s,2) - 
        2*Power(mK,2)*(Power(mKstar,2) + s)))/(2.*Power(mKstar,2)*Power(-Power(mpion,2) + s,2));

  double MaMb = (Power(C,2)*Power(gk,4)*(2*Power(mKstar,6) - 3*Power(mKstar,4)*s + Power(s,3) + 
       Power(mpion,2)*(-Power(mKstar,2) + s)*(-Power(mKstar,2) + s - t) + 
       Power(mK,4)*(Power(mKstar,2) + 4*Power(mpion,2) + s - t) - 3*Power(mKstar,4)*t + Power(mKstar,2)*s*t + 
       2*Power(s,2)*t + Power(mKstar,2)*Power(t,2) + s*Power(t,2) - 
       Power(mK,2)*(3*Power(mKstar,4) + 2*Power(s,2) + Power(mpion,2)*(7*Power(mKstar,2) + 5*s - t) + s*t + 
          Power(t,2) - Power(mKstar,2)*(s + 4*t))))/
   (16.*Power(mKstar,2)*(-Power(mpion,2) + s)*(-Power(mK,2) - Power(mpion,2) + s + t));

  double MaMc = (Power(C,2)*Power(gk,4)*(Power(mKstar,6) - 2*Power(mKstar,4)*s + Power(mKstar,2)*Power(s,2) - 2*Power(mKstar,4)*t + 
       Power(mKstar,2)*s*t + Power(s,2)*t + Power(mKstar,2)*Power(t,2) + s*Power(t,2) - 
       Power(mpion,2)*(-Power(mKstar,2) + s)*(-Power(mKstar,2) + s + t) - 
       Power(mK,2)*(-Power(mKstar,2) - Power(mpion,2) + t)*(-Power(mKstar,2) + s + t)))/
   (8.*Power(mKstar,2)*(-Power(mpion,2) + s)*(-Power(mK,2) + t));

  double MaMd = (-3*Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + Power(mpion,2)*(Power(mKstar,2) - s) + Power(mKstar,2)*s + 
       Power(mK,2)*(Power(mKstar,2) + Power(mpion,2) - t) + Power(mKstar,2)*t + s*t))/
   (16.*Power(mKstar,2)*(-Power(mpion,2) + s));

  double MbMb = -(Power(C,2)*Power(gk,4)*(-Power(mK,6) + 4*Power(mKstar,6) - Power(mpion,6) - 4*Power(mKstar,4)*s - 
        3*Power(mKstar,2)*Power(s,2) + Power(s,3) + Power(mK,4)*(Power(mKstar,2) + 5*Power(mpion,2) + 3*s - t) - 
        4*Power(mKstar,4)*t + 2*Power(mKstar,2)*s*t + 3*Power(s,2)*t - 3*Power(mKstar,2)*Power(t,2) + 
        3*s*Power(t,2) + Power(t,3) + Power(mpion,4)*(Power(mKstar,2) - s + 3*t) + 
        Power(mpion,2)*(-4*Power(mKstar,4) + Power(s,2) + Power(mKstar,2)*(6*s - 2*t) - 2*s*t - 3*Power(t,2)) - 
        Power(mK,2)*(4*Power(mKstar,4) - 5*Power(mpion,4) + 3*Power(s,2) + 2*Power(mKstar,2)*(s - 3*t) + 2*s*t - 
           Power(t,2) + 6*Power(mpion,2)*(Power(mKstar,2) + s + t))))/
   (32.*Power(mKstar,2)*Power(-Power(mK,2) - Power(mpion,2) + s + t,2));

  double MbMc = -(Power(C,2)*Power(gk,4)*(2*Power(mKstar,6) - 3*Power(mKstar,4)*s + Power(mKstar,2)*Power(s,2) - 
        3*Power(mKstar,4)*t + Power(mKstar,2)*s*t + Power(s,2)*t + 2*s*Power(t,2) + Power(t,3) + 
        Power(mpion,4)*(Power(mKstar,2) - s + t) + 
        Power(mK,2)*(4*Power(mpion,4) + Power(mpion,2)*(-7*Power(mKstar,2) + s - 5*t) - 
           (Power(mKstar,2) + s - t)*(-Power(mKstar,2) + t)) - 
        Power(mpion,2)*(3*Power(mKstar,4) + Power(s,2) + s*t + 2*Power(t,2) - Power(mKstar,2)*(4*s + t))))/
   (32.*Power(mKstar,2)*(-Power(mK,2) + t)*(-Power(mK,2) - Power(mpion,2) + s + t));

  double MbMd =(-3*Power(C,2)*Power(gk,4)*(Power(mK,4) - Power(mpion,4) + 4*Power(mKstar,2)*s + Power(s,2) - 
       2*Power(mK,2)*(-2*Power(mKstar,2) + s) - 4*Power(mKstar,2)*t - Power(t,2) + 
       2*Power(mpion,2)*(-2*Power(mKstar,2) + t)))/(64.*Power(mKstar,2)*(-Power(mK,2) - Power(mpion,2) + s + t));

  double McMc = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mpion,4) + Power(-Power(mKstar,2) + t,2) - 
        2*Power(mpion,2)*(Power(mKstar,2) + t)))/(8.*Power(mKstar,2)*Power(-Power(mK,2) + t,2));

  double McMd = (-3*Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + Power(mpion,2)*(Power(mKstar,2) - s) + Power(mKstar,2)*s + 
       Power(mK,2)*(Power(mKstar,2) + Power(mpion,2) - t) + Power(mKstar,2)*t + s*t))/
   (32.*Power(mKstar,2)*(-Power(mK,2) + t));

  double MdMd = (27*Power(C,2)*Power(gk,4))/32.;

  double result = MaMa + 2*(MaMb + MaMc + MaMd) + MbMb + 2*(MbMc + MbMd) + McMc + 2*McMd + MdMd;

  return(result);
}

double C7p2(double s, double t)
{
  double C = 0.059;
  double gk = 9.2;

  double MaMa = -(Power(C,2)*Power(gk,4)*Power(mK,2)*(Power(mpion,4) + Power(-Power(mKstar,2) + t,2) - 
        2*Power(mpion,2)*(Power(mKstar,2) + t)))/(16.*Power(mKstar,2)*Power(-Power(mK,2) + t,2));

  double MaMb = (Power(C,2)*Power(gk,4)*(2*Power(mKstar,6) - 3*Power(mKstar,4)*s + Power(mKstar,2)*Power(s,2) - 
       3*Power(mKstar,4)*t + Power(mKstar,2)*s*t + Power(s,2)*t + 2*s*Power(t,2) + Power(t,3) + 
       Power(mpion,4)*(Power(mKstar,2) - s + t) + 
       Power(mK,2)*(4*Power(mpion,4) + Power(mpion,2)*(-7*Power(mKstar,2) + s - 5*t) - 
          (Power(mKstar,2) + s - t)*(-Power(mKstar,2) + t)) - 
       Power(mpion,2)*(3*Power(mKstar,4) + Power(s,2) + s*t + 2*Power(t,2) - Power(mKstar,2)*(4*s + t))))/
   (64.*Power(mKstar,2)*(-Power(mK,2) + t)*(-Power(mK,2) - Power(mpion,2) + s + t));

  double MaMc = -(Power(C,2)*Power(gk,4)*(-Power(mKstar,4) + Power(mpion,2)*(Power(mKstar,2) - s) + Power(mKstar,2)*s + 
        Power(mK,2)*(Power(mKstar,2) + Power(mpion,2) - t) + Power(mKstar,2)*t + s*t))/
   (64.*Power(mKstar,2)*(-Power(mK,2) + t));

  double MbMb = -(Power(C,2)*Power(gk,4)*(-Power(mK,6) + 4*Power(mKstar,6) - Power(mpion,6) - 4*Power(mKstar,4)*s - 
        3*Power(mKstar,2)*Power(s,2) + Power(s,3) + Power(mK,4)*(Power(mKstar,2) + 5*Power(mpion,2) + 3*s - t) - 
        4*Power(mKstar,4)*t + 2*Power(mKstar,2)*s*t + 3*Power(s,2)*t - 3*Power(mKstar,2)*Power(t,2) + 
        3*s*Power(t,2) + Power(t,3) + Power(mpion,4)*(Power(mKstar,2) - s + 3*t) + 
        Power(mpion,2)*(-4*Power(mKstar,4) + Power(s,2) + Power(mKstar,2)*(6*s - 2*t) - 2*s*t - 3*Power(t,2)) - 
        Power(mK,2)*(4*Power(mKstar,4) - 5*Power(mpion,4) + 3*Power(s,2) + 2*Power(mKstar,2)*(s - 3*t) + 2*s*t - 
           Power(t,2) + 6*Power(mpion,2)*(Power(mKstar,2) + s + t))))/
   (64.*Power(mKstar,2)*Power(-Power(mK,2) - Power(mpion,2) + s + t,2));
  
  double MbMc = (Power(C,2)*Power(gk,4)*(Power(mK,4) - Power(mpion,4) + 4*Power(mKstar,2)*s + Power(s,2) - 
       2*Power(mK,2)*(-2*Power(mKstar,2) + s) - 4*Power(mKstar,2)*t - Power(t,2) + 
       2*Power(mpion,2)*(-2*Power(mKstar,2) + t)))/(128.*Power(mKstar,2)*(-Power(mK,2) - Power(mpion,2) + s + t));

  double McMc = (3*Power(C,2)*Power(gk,4))/64.;

  double result = MaMa + 2*(MaMb + MaMc) + MbMb + 2*MbMc + McMc;
  return(result);
}

inline double Power(double x, int a)
{
    double result = 1.0;
    for(int i=0; i<a; i++)
       result *= x;
    return(result);
}
