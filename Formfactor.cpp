#include <cmath>
#include <iostream>

#include "parameters.h"
#include "Formfactor.h"
using namespace std;

void Calculate_Formfactor(double* E, double* Formfactor, int npoint, int channel)
{
   double tbar;
   const double Lambda = 1.0;
   switch(channel)
   {
      case 1: //C1: pi + rho -> pi + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = 34.5096*pow(E[i], 0.737) - 67.557*pow(E[i], 0.7584) + 32.858*pow(E[i], 0.7806); 
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      case 2: //C1_omega: pi + rho -> omega -> pi + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = -61.595*pow(E[i], 0.9979) + 28.592*pow(E[i], 1.1579) + 37.738*pow(E[i], 0.9317) - 5.282*pow(E[i], 1.3686);
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      case 3: //C2: pi + pi -> rho + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = 34.5096*pow(E[i], 0.737) - 67.557*pow(E[i], 0.7584) + 32.858*pow(E[i], 0.7806); 
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      case 4: //C4: pi + Kstar -> K + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = 34.5096*pow(E[i], 0.737) - 67.557*pow(E[i], 0.7584) + 32.858*pow(E[i], 0.7806); 
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      case 5: //C5: pi + K -> Kstar + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = 34.5096*pow(E[i], 0.737) - 67.557*pow(E[i], 0.7584) + 32.858*pow(E[i], 0.7806); 
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      case 6: //C6: rho + K -> K + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = -76.424*pow(E[i], 0.6236) + 36.944*pow(E[i], 0.6604) + 39.0448*pow(E[i], 0.5873);
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      case 7: //C7: K + Kstar -> pi + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = -76.424*pow(E[i], 0.6236) + 36.944*pow(E[i], 0.6604) + 39.0448*pow(E[i], 0.5873);
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      case 8: //C3: rho -> pi + pi + gamma
         for(int i=0; i<npoint; i++)
         {
            tbar = 34.5096*pow(E[i], 0.737) - 67.557*pow(E[i], 0.7584) + 32.858*pow(E[i], 0.7806); 
            Formfactor[i] = pow((2*Lambda*Lambda/(2*Lambda*Lambda - tbar)), 8); 
         }
         break;
      default:
         cout << "calcualte Formfactor ERROR: can not find the corresponding channel! (channel =  " << channel << ")" << endl;
         exit(1);
   }
}
