#ifndef HG_1to3_decay_H
#define HG_1to3_decay_H

#include <string>
#include "ParameterReader.h"
#include "Physicalconstants.h"
#include "chemical_potential.h"

class HG_1to3_decay;

struct CCallbackHolder
{
   HG_1to3_decay* clsPtr;
   void *params;
};

class HG_1to3_decay
{
   private:
      ParameterReader *paraRdr;

      double eps;
      int n_Eq, n_Temp;
      double *Eq_tb, *T_tb;
      double **equilibrium_results, **viscous_results, **bulkvis_results;

      int channel;
      double *m, *mu;
      string filename;
      
      // Gaussian quadrature points for phase space integrations 
      int n_s;
      double *s_pt, *s_weight, *s_pt_standard, *s_weight_standard;
      int n_t;
      double **t_pt, **t_weight, **t_pt_standard, **t_weight_standard;
      double **Matrix_elements_sq_ptr;
      int n_E1;
      double *E1_pt_standard, *E1_weight_standard;
      int n_E2;
      double **E2_pt_standard, **E2_weight_standard;

      double deltaf_alpha;

   public:
      HG_1to3_decay(ParameterReader* paraRdr_in);
      ~HG_1to3_decay();
      
      void output_emissionrateTable();
      int Calculate_emissionrates(Chemical_potential* chempotential_ptr, int channel, string filename_in);
      void set_particleMass();
      void set_gausspoints();
      double Integrate_E1(double Eq, double T, double s, double t, double* results);
      double Integrate_E2(double Eq, double T, double s, double t, double E1, double* results);
      double viscous_integrand(double s, double t, double E1, double E2, double Eq, double T, double f0_E1, double f0_E2, double f0_E3);
      double bulkvis_integrand(double s, double t, double E1, double E2, double Eq, double T, double f0_E1, double f0_E2, double f0_E3);
      double Bose_distribution(double E, double T, double mu);
      double deltaf_chi(double p);
      
      //static call back function for gsl integration
      double Rateintegrands(double s, void *params);
      static double CCallback_Rateintegrands(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->Rateintegrands(x, h->params);
      }
      double Rateintegrandt(double t, void *params);
      static double CCallback_Rateintegrandt(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->Rateintegrandt(x, h->params);
      }
      double RateintegrandE1(double E1, void *params);
      static double CCallback_RateintegrandE1(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->RateintegrandE1(x, h->params);
      }
      double RateintegrandE2(double E2, void *params);
      static double CCallback_RateintegrandE2(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->RateintegrandE2(x, h->params);
      }

      //matrix elements squared
      double Matrix_elements_sq(double s, double t);
      double Matrix_elements_sq_C1(double s, double t);
      double Matrix_elements_sq_C1_omega(double s, double t);
      double Matrix_elements_sq_C2(double s, double t);
      double Matrix_elements_sq_C3(double s, double t);
      double Matrix_elements_sq_C4(double s, double t);
      double Matrix_elements_sq_C5(double s, double t);
      double Matrix_elements_sq_C6(double s, double t);
      double Matrix_elements_sq_C7(double s, double t);
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
      
      double Power(double x, int a);
};

#endif
