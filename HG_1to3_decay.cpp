#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_integration.h>

#include "HG_1to3_decay.h"
#include "ParameterReader.h"
#include "Physicalconstants.h"
#include "gauss_quadrature.h"
#include "Arsenal.h"
#include "Formfactor.h"

using namespace std;

HG_1to3_decay::HG_1to3_decay(ParameterReader* paraRdr_in)
{
   eps = 1e-16;
   paraRdr = paraRdr_in;

   n_Eq = paraRdr->getVal("n_Eq");
   double Eq_i = paraRdr->getVal("Eq_min");
   double dEq = paraRdr->getVal("dEq");
   Eq_tb = new double [n_Eq];
   for(int i=0; i<n_Eq; i++)
      Eq_tb[i] = Eq_i + i*dEq;

   n_Temp = paraRdr->getVal("n_Temp");
   double T_i = paraRdr->getVal("T_min");
   double dT = paraRdr->getVal("dT");
   T_tb = new double [n_Temp];
   for(int i=0; i<n_Temp; i++)
      T_tb[i] = T_i + i*dT;

   equilibrium_results = new double* [n_Eq];
   viscous_results = new double*[n_Eq];
   bulkvis_results = new double*[n_Eq];
   for(int i=0; i<n_Eq; i++)
   {
      equilibrium_results[i] = new double [n_Temp];
      viscous_results[i] = new double [n_Temp];
      bulkvis_results[i] = new double [n_Temp];
   }
   
   //initialize the Gaussian quadrature lattices
   n_s = paraRdr->getVal("n_s");
   s_pt = new double [n_s];
   s_weight = new double [n_s];
   s_pt_standard = new double [n_s];
   s_weight_standard = new double [n_s];
   
   t_pt = new double* [n_s];
   t_weight = new double* [n_s];
   t_pt_standard = new double* [n_s];
   t_weight_standard = new double* [n_s];
   Matrix_elements_sq_ptr = new double* [n_s];
   n_t = paraRdr->getVal("n_t");
   for(int i=0; i<n_s; i++)
   {
       t_pt[i] = new double [n_t];
       t_weight[i] = new double [n_t];
       t_pt_standard[i] = new double [n_t];
       t_weight_standard[i] = new double [n_t];
       Matrix_elements_sq_ptr[i] = new double [n_t];
   }
   
   n_E1 = paraRdr->getVal("n_E1");
   E1_pt_standard = new double [n_E1];
   E1_weight_standard = new double [n_E1];
   n_E2 = paraRdr->getVal("n_E2");
   E2_pt_standard = new double* [3];
   E2_weight_standard = new double* [3];
   for(int i=0; i<3; i++)
   {
      E2_pt_standard[i] = new double [n_E2];
      E2_weight_standard[i] = new double [n_E2];
   }
   
   m = new double [3];
   mu = new double [3];
   deltaf_alpha = paraRdr->getVal("deltaf_alpha");
   bulk_deltaf_kind = paraRdr->getVal("bulk_deltaf_kind");

   if(bulk_deltaf_kind == 0)
       bulkdf_coeff = new Table ("chemical_potential_tb/s95p/s95p-PCE165-v0/BulkDf_Coefficients_Hadrons_CE.dat");
}

HG_1to3_decay::~HG_1to3_decay()
{
   delete[] Eq_tb;
   delete[] T_tb;
   for(int i=0; i<n_Eq; i++)
   {
      delete[] equilibrium_results[i];
      delete[] viscous_results[i];
      delete[] bulkvis_results[i];
   }
   delete[] equilibrium_results;
   delete[] viscous_results;
   delete[] bulkvis_results;

   delete[] s_pt;
   delete[] s_weight;
   delete[] s_pt_standard;
   delete[] s_weight_standard;
   for(int i=0; i<n_s; i++)
   {
      delete[] t_pt[i];
      delete[] t_weight[i];
      delete[] t_pt_standard[i];
      delete[] t_weight_standard[i];
      delete[] Matrix_elements_sq_ptr[i];
   }
   delete[] t_pt;
   delete[] t_weight;
   delete[] t_pt_standard;
   delete[] t_weight_standard;
   delete[] Matrix_elements_sq_ptr;
   
   delete[] E1_pt_standard;
   delete[] E1_weight_standard;
   for(int i=0; i<3; i++)
   {
      delete[] E2_pt_standard[i];
      delete[] E2_weight_standard[i];
   }
   delete[] E2_pt_standard;
   delete[] E2_weight_standard;
   delete[] m;
   delete[] mu;

   if(bulk_deltaf_kind == 0)
       delete bulkdf_coeff;

}

void HG_1to3_decay::output_emissionrateTable()
{
   ostringstream output_file_eqrate;
   ostringstream output_file_viscous;
   ostringstream output_file_bulkvis;
   output_file_eqrate << "rate_" << filename << "_eqrate.dat";
   output_file_viscous << "rate_" << filename << "_viscous.dat";
   output_file_bulkvis << "rate_" << filename << "_bulkvis.dat";
   ofstream of_eqrate(output_file_eqrate.str().c_str());
   ofstream of_viscous(output_file_viscous.str().c_str());
   ofstream of_bulkvis(output_file_bulkvis.str().c_str());
   for(int j=0; j<n_Temp; j++)
   {
      for(int i=0; i<n_Eq; i++)
      {
         of_eqrate << scientific << setw(20) << setprecision(8)
                   << equilibrium_results[i][j] << "   ";
         of_viscous << scientific << setw(20) << setprecision(8)
                    << viscous_results[i][j] << "   ";
         of_bulkvis << scientific << setw(20) << setprecision(8)
                    << bulkvis_results[i][j] << "   ";
      }
      of_eqrate << endl;
      of_viscous << endl;
      of_bulkvis << endl;
   }
}

int HG_1to3_decay::Calculate_emissionrates(Chemical_potential* chempotential_ptr, int channel_in, string filename_in)
{
   double* results = new double [3];

   filename = filename_in; 
   channel = channel_in;

   set_particleMass();
   set_gausspoints();
   
   // calculate matrix elements squared
   for(int i=0; i<n_s; i++)
      for(int j=0; j<n_t; j++)
           Matrix_elements_sq_ptr[i][j] = Matrix_elements_sq(s_pt[i], t_pt[i][j]);

   double* Formfactor_tb = new double [n_Eq];
   double* mu1_tb = new double [n_Temp];
   double* mu2_tb = new double [n_Temp];
   double* mu3_tb = new double [n_Temp];

   // calculate form factor
   Calculate_Formfactor(Eq_tb, Formfactor_tb, n_Eq, channel);
   // calculate chemical potentials
   chempotential_ptr->Calculate_mu(T_tb, mu1_tb, mu2_tb, mu3_tb, n_Temp, channel);

   double Eq;
   double formfactor;
   double T;
   for(int j=0; j<n_Temp; j++)
   {
      T = T_tb[j];
      mu[0] = mu1_tb[j];
      mu[1] = mu2_tb[j];
      mu[2] = mu3_tb[j];
      for(int i=0; i<n_Eq; i++)
      {
          Eq = Eq_tb[i];
          formfactor = Formfactor_tb[i];
          double prefactor = 1./16./pow(2.0*M_PI, 7)/Eq*formfactor;
          
          double equilibrium_result_s = 0.0;
          double viscous_result_s = 0.0;
          double bulkvis_result_s = 0.0;
          for(int k=0; k<n_s; k++)
          {
             double equilibrium_result_t = 0.0;
             double viscous_result_t = 0.0;
             double bulkvis_result_t = 0.0;
             for(int l=0; l<n_t; l++)
             {
                Integrate_E1(Eq, T, s_pt[k], t_pt[k][l], results);
                equilibrium_result_t += Matrix_elements_sq_ptr[k][l]*results[0]*t_weight[k][l];
                viscous_result_t += Matrix_elements_sq_ptr[k][l]*results[1]*t_weight[k][l];
                bulkvis_result_t += Matrix_elements_sq_ptr[k][l]*results[2]*t_weight[k][l];
             }
             equilibrium_result_s += equilibrium_result_t*s_weight[k];
             viscous_result_s += viscous_result_t*s_weight[k];
             bulkvis_result_s += bulkvis_result_t*s_weight[k];
          }
          
          equilibrium_results[i][j] = equilibrium_result_s*prefactor/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
          viscous_results[i][j] = viscous_result_s*prefactor/(Eq*Eq)/pow(hbarC, 4); // convert units to 1/(GeV^4 fm^4) for the emission rates
          bulkvis_results[i][j] = bulkvis_result_s*prefactor/pow(hbarC, 4); // convert units to 1/(GeV^4 fm^4) for the emission rates
      }
   }
   output_emissionrateTable();

   delete[] mu1_tb;
   delete[] mu2_tb;
   delete[] mu3_tb;
   delete[] Formfactor_tb;
   delete [] results;
   return 0;
}


void HG_1to3_decay::set_gausspoints()
{
   if(m[0] < (m[1] + m[2]))
   {
      cout << "Error: decay particle mass is smaller than daughter particles! Please check!" << endl;
      exit(0);
   }
   double s_min = m[2]*m[2];
   double s_max = (m[0]-m[1])*(m[0]-m[1]);
  
   gauss_quadrature(n_s, 1, 0.0, 0.0, s_min, s_max, s_pt, s_weight);
  
   for(int i=0; i<n_s; i++)
   {
      double s = s_pt[i];
      double t_min;
      double t_max;
      t_min = m[0]*m[0] + m[2]*m[2] - 2*(s + m[0]*m[0] - m[1]*m[1])*(s+m[2]*m[2])
              /4/s - 2*sqrt((s+m[0]*m[0]-m[1]*m[1])*(s+m[0]*m[0]-m[1]*m[1]) 
              - 4*s*m[0]*m[0])*(s - m[2]*m[2])/4/s;
      t_max = m[0]*m[0] + m[2]*m[2] - 2*(s + m[0]*m[0] - m[1]*m[1])*(s+m[2]*m[2])
              /4/s + 2*sqrt((s+m[0]*m[0]-m[1]*m[1])*(s+m[0]*m[0]-m[1]*m[1])
              - 4*s*m[0]*m[0])*(s - m[2]*m[2])/4/s;

      gauss_quadrature(n_t, 1, 0.0, 0.0, t_min, t_max, t_pt[i], t_weight[i]);
    }
    
    gauss_quadrature_standard(n_E1, 5, 0.0, 0.0, 0.0, 1.0, E1_pt_standard, E1_weight_standard);

    // use Chebyshev–Gauss quadrature for channels: pi + rho, pi + Kstar, rho + K, and K + Kstar
    gauss_quadrature_standard(n_E2, 2, 0.0, 0.0, 0.0, 1.0, E2_pt_standard[0], E2_weight_standard[0]);

    // use Jacobi-Gauss quadrature for channels: pi + pi, pi + K
    gauss_quadrature_standard(n_E2, 4, 0.0, -0.5, 0.0, 1.0, E2_pt_standard[1], E2_weight_standard[1]);
    gauss_quadrature_standard(n_E2, 4, -0.5, 0.0, 0.0, 1.0, E2_pt_standard[2], E2_weight_standard[2]);

    return;
}

void HG_1to3_decay::set_particleMass()
{
   if(channel == 1)
   {
      m[0] = mrho;
      m[1] = mpion;
      m[2] = mpion;
   }
   else if (channel == 2)
   {
      m[0] = mrho;
      m[1] = mpion;
      m[2] = mpion;
   }
   else if (channel == 3)
   {
      m[0] = mpion;
      m[1] = mpion;
      m[2] = mrho;
   }
   else if (channel == 4)
   {
      m[0] = mKstar;
      m[1] = mpion;
      m[2] = mK;
   }
   else if (channel == 5) 
   {
      m[0] = mK;
      m[1] = mpion;
      m[2] = mKstar;
   }
   else if (channel == 6)
   {
      m[0] = mrho;
      m[1] = mK;
      m[2] = mK;
   }
   else if (channel == 7)
   {
      m[0] = mKstar;
      m[1] = mK;
      m[2] = mpion;
   }
   else if (channel == 8)
   {
      m[0] = mrho;
      m[1] = mpion;
      m[2] = mpion;
   }
   else
   {
      cout << "Error:: set_particleMass: input channel is invalid, channel = " << channel << endl;
      exit(1);
   }
   return;
}


double HG_1to3_decay::Integrate_E1(double Eq, double T, double s, double t, double* results)
{
   double equilibrium_result = 0.0e0;
   double viscous_result = 0.0e0;
   double bulkvis_result = 0.0e0;
   double E1_min;
   double u = - s - t + m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
   E1_min = Eq*m[0]*m[0]/(m[0]*m[0] - u) + (m[0]*m[0] - u)/4/Eq;

   double* E1_pt = new double [n_E1];
   double* E1_weight = new double [n_E1];
   for(int i=0; i<n_E1; i++)
   {
      E1_pt[i] = E1_pt_standard[i];
      E1_weight[i] = E1_weight_standard[i];
   }
   
   double slope = 1./T;
   scale_gausspoints(n_E1, 5, 0.0, 0.0, E1_min, slope, E1_pt, E1_weight);

   for(int i=0; i<n_E1; i++)
   {
      Integrate_E2(Eq, T, s, t, E1_pt[i], results);
      equilibrium_result += results[0]*E1_weight[i];
      viscous_result += results[1]*E1_weight[i];
      bulkvis_result += results[2]*E1_weight[i];
   }
   
   results[0] = equilibrium_result;
   results[1] = viscous_result;
   results[2] = bulkvis_result;

   delete[] E1_pt;
   delete[] E1_weight;

   return(0);
}

double HG_1to3_decay::Integrate_E2(double Eq, double T, double s, double t, double E1, double* results)
{
   double equilibrium_result = 0.0;
   double viscous_result = 0.0;
   double bulkvis_result = 0.0;
   double E2_min;
   double E2_max;
   double min_1 = Eq*m[1]*m[1]/(t - m[1]*m[1]) + (t - m[1]*m[1])/4/Eq;

   double a = - (s + t - m[1]*m[1] - m[2]*m[2])*(s + t - m[1]*m[1] - m[2]*m[2]);
   double b = - Eq*((s + t - m[1]*m[1] - m[2]*m[2])*(s - m[0]*m[0] - m[1]*m[1]) 
              - 2*m[0]*m[0]*(m[1]*m[1] - t)) - E1*(m[1]*m[1] - t)
              *(s + t - m[1]*m[1] - m[2]*m[2]);
   double c = - (t - m[1]*m[1])*(t - m[1]*m[1])*E1*E1
              - 2*Eq*(2*m[1]*m[1]*(s + t - m[1]*m[1] - m[2]*m[2]) 
              - (m[1]*m[1] - t)*(s - m[0]*m[0] - m[1]*m[1]))*E1
              + 4*Eq*Eq*m[0]*m[0]*m[1]*m[1] + m[1]*m[1]*(s + t - m[1]*m[1] 
              - m[2]*m[2])*(s + t - m[1]*m[1] - m[2]*m[2]) + m[0]*m[0]
              *(m[1]*m[1] -t)*(m[1]*m[1] -t)
              - Eq*Eq*(s - m[0]*m[0] - m[1]*m[1])*(s - m[0]*m[0] - m[1]*m[1])
              + (s - m[0]*m[0] - m[1]*m[1])*(t - m[1]*m[1])*(s + t - m[1]*m[1] 
              - m[2]*m[2]);

   if((b*b - a*c) >= 0) 
   {
      double min_2 = (-b + sqrt(b*b - a*c))/a;
      if(min_1 < min_2)
         E2_min = min_2;
      else
         E2_min = min_1;
      E2_max = (-b - sqrt(b*b - a*c))/a;

      if(E2_max < E2_min)
      {
         results[0] = 0.0e0;
         results[1] = 0.0e0;
         results[2] = 0.0e0;
         return (0.0);
      }
   
      double mu1 = mu[0];
      double mu2 = mu[1];
      double mu3 = mu[2];
      double common_factor;

      double* E2_pt = new double [n_E2];
      double* E2_weight = new double [n_E2];

      double* bulkvis_B0 = new double [3];
      double* bulkvis_D0 = new double [3];
      double* bulkvis_E0 = new double [3];

      if(channel == 3 || channel == 5)
      {
         double E2_cut = E2_min + (E2_max - E2_min)/100.;
         for(int i=0; i<n_E2; i++)
         {
            E2_pt[i] = E2_pt_standard[1][i];
            E2_weight[i] = E2_weight_standard[1][i];
         }
         scale_gausspoints(n_E2, 4, 0.0, -0.5, E2_min, E2_cut, E2_pt, E2_weight);
         for(int i=0; i<n_E2; i++)
         {
            double f0_E1 = Bose_distribution(E1, T, mu1);
            double f0_E2 = Bose_distribution(E2_pt[i], T, mu2);
            double f0_E3 = Bose_distribution(E1 + E2_pt[i] - Eq, T, mu3);
            common_factor = f0_E1*(1 + f0_E2)*(1 + f0_E3)/(sqrt(a*E2_pt[i]*E2_pt[i] + 2*b*E2_pt[i] + c));
            equilibrium_result += common_factor*1.*E2_weight[i];
            viscous_result += common_factor*viscous_integrand(s, t, E1, E2_pt[i], Eq, T, f0_E1, f0_E2, f0_E3)*E2_weight[i];

            if(bulk_deltaf_kind == 0)
            {
                get_bulkvis_coefficients_14moment(T, bulkvis_B0, bulkvis_D0, bulkvis_E0);
                bulkvis_result += common_factor*bulkvis_integrand_14moment(E1, E2_pt[i], Eq, f0_E1, f0_E2, f0_E3, bulkvis_B0, bulkvis_D0, bulkvis_E0)*E2_weight[i];
            }
            else if (bulk_deltaf_kind == 1)
            {
                double bulkvis_Cbulk, bulkvis_e2;
                get_bulkvis_coefficients_relaxation(T, &bulkvis_Cbulk, &bulkvis_e2);
                bulkvis_result += common_factor*bulkvis_integrand_relaxation(T, E1, E2_pt[i], Eq, f0_E1, f0_E2, f0_E3, bulkvis_Cbulk, bulkvis_e2)*E2_weight[i];
            }
         }

         for(int i=0; i<n_E2; i++)
         {
            E2_pt[i] = E2_pt_standard[2][i];
            E2_weight[i] = E2_weight_standard[2][i];
         }
         scale_gausspoints(n_E2, 4, -0.5, 0.0, E2_cut, E2_max, E2_pt, E2_weight);
         for(int i=0; i<n_E2; i++)
         {
            double f0_E1 = Bose_distribution(E1, T, mu1);
            double f0_E2 = Bose_distribution(E2_pt[i], T, mu2);
            double f0_E3 = Bose_distribution(E1 + E2_pt[i] - Eq, T, mu3);
            common_factor = f0_E1*(1 + f0_E2)*(1 + f0_E3)/(sqrt(a*E2_pt[i]*E2_pt[i] + 2*b*E2_pt[i] + c));
            equilibrium_result += common_factor*1.*E2_weight[i];
            viscous_result += common_factor*viscous_integrand(s, t, E1, E2_pt[i], Eq, T, f0_E1, f0_E2, f0_E3)*E2_weight[i];
            if(bulk_deltaf_kind == 0)
            {
                get_bulkvis_coefficients_14moment(T, bulkvis_B0, bulkvis_D0, bulkvis_E0);
                bulkvis_result += common_factor*bulkvis_integrand_14moment(E1, E2_pt[i], Eq, f0_E1, f0_E2, f0_E3, bulkvis_B0, bulkvis_D0, bulkvis_E0)*E2_weight[i];
            }
            else if(bulk_deltaf_kind == 1)
            {
                double bulkvis_Cbulk, bulkvis_e2;
                get_bulkvis_coefficients_relaxation(T, &bulkvis_Cbulk, &bulkvis_e2);
                bulkvis_result += common_factor*bulkvis_integrand_relaxation(T, E1, E2_pt[i], Eq, f0_E1, f0_E2, f0_E3, bulkvis_Cbulk, bulkvis_e2)*E2_weight[i];
            }
         }
      }
      else
      {
         for(int i=0; i<n_E2; i++)
         {
            E2_pt[i] = E2_pt_standard[0][i];
            E2_weight[i] = E2_weight_standard[0][i];
         }
         scale_gausspoints(n_E2, 2, 0.0, 0.0, E2_min, E2_max, E2_pt, E2_weight);
         for(int i=0; i<n_E2; i++)
         {
            double f0_E1 = Bose_distribution(E1, T, mu1);
            double f0_E2 = Bose_distribution(E2_pt[i], T, mu2);
            double f0_E3 = Bose_distribution(E1 - E2_pt[i] - Eq, T, mu3);
            common_factor = f0_E1*(1 + f0_E2)*(1 + f0_E3)/(sqrt(a*E2_pt[i]*E2_pt[i] + 2*b*E2_pt[i] + c));
            equilibrium_result += common_factor*1.*E2_weight[i];
            viscous_result += common_factor*viscous_integrand(s, t, E1, E2_pt[i], Eq, T, f0_E1, f0_E2, f0_E3)*E2_weight[i];
            if(bulk_deltaf_kind == 0)
            {
                get_bulkvis_coefficients_14moment(T, bulkvis_B0, bulkvis_D0, bulkvis_E0);
                bulkvis_result += common_factor*bulkvis_integrand_14moment(E1, E2_pt[i], Eq, f0_E1, f0_E2, f0_E3, bulkvis_B0, bulkvis_D0, bulkvis_E0)*E2_weight[i];
            }
            else if(bulk_deltaf_kind == 1)
            {
                double bulkvis_Cbulk, bulkvis_e2;
                get_bulkvis_coefficients_relaxation(T, &bulkvis_Cbulk, &bulkvis_e2);
                bulkvis_result += common_factor*bulkvis_integrand_relaxation(T, E1, E2_pt[i], Eq, f0_E1, f0_E2, f0_E3, bulkvis_Cbulk, bulkvis_e2)*E2_weight[i];
            }
         }
      }

      delete[] E2_pt;
      delete[] E2_weight;
      delete[] bulkvis_B0;
      delete[] bulkvis_D0;
      delete[] bulkvis_E0;
   }
   else  // no kinematic phase space
   {
      equilibrium_result = 0.0e0;
      viscous_result = 0.0e0;
      bulkvis_result = 0.0e0;
   }
   results[0] = equilibrium_result;
   results[1] = viscous_result;
   results[2] = bulkvis_result;

   return(0);
}

double HG_1to3_decay::viscous_integrand(double s, double t, double E1, double E2, double Eq, double T, double f0_E1, double f0_E2, double f0_E3)
{
   double m1 = m[0];
   double m2 = m[1];
   double m3 = m[2];
   double E3 = E1 - E2 - Eq;
   double p1 = sqrt(E1*E1 - m1*m1);
   double p2 = sqrt(E2*E2 - m2*m2);
   double p3 = sqrt(E3*E3 - m3*m3);
   double costheta1 = (- s - t + m2*m2 + m3*m3 + 2*E1*Eq)/(2*p1*Eq);
   double costheta2 = ( - t + m2*m2 + 2*E2*Eq)/(2*p2*Eq);
   double p3_z = p1*costheta1 - p2*costheta2 - Eq; 

   double integrand =   (1. + f0_E1)*deltaf_chi(p1/T)*0.5*(-1. + 3.*costheta1*costheta1) 
                      + f0_E2*deltaf_chi(p2/T)*0.5*(-1. + 3.*costheta2*costheta2) 
                      + f0_E3*deltaf_chi(p3/T)/p3/p3*(-0.5*p3*p3 + 1.5*p3_z*p3_z);

   return(integrand);
}

void HG_1to3_decay::get_bulkvis_coefficients_14moment(double T, double* bulkvis_B0, double* bulkvis_D0, double * bulkvis_E0)
{
   double T_fm = T/hbarC;  // convert to [1/fm]

   for(int i = 0; i < 3; i++)
   {
      bulkvis_B0[i] = bulkdf_coeff->interp(1, 2, T_fm, 5)/pow(hbarC, 3);  // [fm^3/GeV^3]
      bulkvis_D0[i] = bulkdf_coeff->interp(1, 3, T_fm, 5)/pow(hbarC, 2);  // [fm^3/GeV^2]
      bulkvis_E0[i] = bulkdf_coeff->interp(1, 4, T_fm, 5)/pow(hbarC, 3);  // [fm^3/GeV^3]

      // parameterization for mu = 0
      //bulkvis_B0[i] = exp(-15.04512474*T_fm + 11.76194266)/pow(hbarC, 3);   // convert to [fm^3/GeV^3]
      //bulkvis_D0[i] = exp( -12.45699277*T_fm + 11.4949293)/pow(hbarC, 2);   // convert to [fm^3/GeV^2]
      //bulkvis_E0[i] = -exp(-14.45087586*T_fm + 11.62716548)/pow(hbarC, 3);  // convert to [fm^3/GeV^3]
   }
   return;
}

void HG_1to3_decay::get_bulkvis_coefficients_relaxation(double T, double* bulkvis_Cbulk, double* bulkvis_e2)
// coefficients for bulk viscous corrections (fits from Gabriel Denicol, derived from relaxation time approximation)
{
   double T_fm = T/hbarC;  // convert to [1/fm]
   double T_power[11];
   T_power[0] = 1.0;
   for(int i = 1; i < 11; i++)
       T_power[i] = T_power[i-1]*T_fm;

   *bulkvis_Cbulk = (  642096.624265727 - 8163329.49562861*T_power[1] 
                     + 47162768.4292073*T_power[2] - 162590040.002683*T_power[3] 
                     + 369637951.096896*T_power[4] - 578181331.809836*T_power[5] 
                     + 629434830.225675*T_power[6] - 470493661.096657*T_power[7] 
                     + 230936465.421*T_power[8] - 67175218.4629078*T_power[9]
                     + 8789472.32652964*T_power[10]);
   *bulkvis_e2 = (  1.18171174036192 - 17.6740645873717*T_power[1]
                  + 136.298469057177*T_power[2] - 635.999435106846*T_power[3]
                  + 1918.77100633321*T_power[4] - 3836.32258307711*T_power[5]
                  + 5136.35746882372*T_power[6] - 4566.22991441914*T_power[7]
                  + 2593.45375240886*T_power[8] - 853.908199724349*T_power[9]
                  + 124.260460450113*T_power[10]);
   return;
}

double HG_1to3_decay::bulkvis_integrand_14moment(double E1, double E2, double Eq, double f0_E1, double f0_E2, double f0_E3, double* bulkvis_B0, double* bulkvis_D0, double* bulkvis_E0)
{
   double E3 = E1 - E2 - Eq;
   double integrand =   (1. + f0_E1)*(-bulkvis_B0[0]*m[0]*m[0] - E1*bulkvis_D0[0] - E1*E1*bulkvis_E0[0])
                      + f0_E2*(-bulkvis_B0[1]*m[1]*m[1] - E2*bulkvis_D0[1] - E2*E2*bulkvis_E0[1])
                      + f0_E3*(-bulkvis_B0[2]*m[2]*m[2] - E3*bulkvis_D0[2] - E3*E3*bulkvis_E0[2]);

   return(integrand);
}

double HG_1to3_decay::bulkvis_integrand_relaxation(double T, double E1, double E2, double Eq, double f0_E1, double f0_E2, double f0_E3, double bulkvis_Cbulk, double bulkvis_e2)
{
   double E3 = E1 - E2 - Eq;
   double E1_over_T = E1/T;
   double E2_over_T = E2/T;
   double E3_over_T = E3/T;
   double integrand =   (1. + f0_E1)*bulkvis_Cbulk/(E1_over_T)*(-m[0]*m[0]/(3.*T*T) + bulkvis_e2*E1_over_T*E1_over_T)
                      + f0_E2*bulkvis_Cbulk/(E2_over_T)*(-m[1]*m[1]/(3.*T*T) + bulkvis_e2*E2_over_T*E2_over_T)
                      + f0_E3*bulkvis_Cbulk/(E3_over_T)*(-m[2]*m[2]/(3.*T*T) + bulkvis_e2*E3_over_T*E3_over_T);
   return(integrand);
}

double HG_1to3_decay::Bose_distribution(double E, double T, double mu)
{
   return(1.0/(exp((E-mu)/T)-1.0));
}

double HG_1to3_decay::deltaf_chi(double p)
{ 
    return(pow(p, deltaf_alpha));
}

double HG_1to3_decay::Rateintegrands(double s, void *params)
{
    double *par = (double*)params;
    double rateType = par[0];
    double Temp = par[1];
    double Eq = par[2];

    double t_min, t_max;
    t_min = m[0]*m[0] + m[2]*m[2] - 2*(s + m[0]*m[0] - m[1]*m[1])*(s+m[2]*m[2])
            /4/s - 2*sqrt((s+m[0]*m[0]-m[1]*m[1])*(s+m[0]*m[0]-m[1]*m[1]) 
            - 4*s*m[0]*m[0])*(s - m[2]*m[2])/4/s;
    t_max = m[0]*m[0] + m[2]*m[2] - 2*(s + m[0]*m[0] - m[1]*m[1])*(s+m[2]*m[2])
            /4/s + 2*sqrt((s+m[0]*m[0]-m[1]*m[1])*(s+m[0]*m[0]-m[1]*m[1])
            - 4*s*m[0]*m[0])*(s - m[2]*m[2])/4/s;

    double *paramsPtr = new double [4];
    paramsPtr[0] = rateType;
    paramsPtr[1] = Temp;
    paramsPtr[2] = Eq;
    paramsPtr[3] = s;
    CCallbackHolder *Callback_params = new CCallbackHolder;
    Callback_params->clsPtr = this;
    Callback_params->params = paramsPtr;
    int maxInteration = 1000;
    gsl_integration_workspace *gsl_workSpace = gsl_integration_workspace_alloc(maxInteration);
    double gslresult, gslerror;
    int status;
    gsl_function gslFunc;
    gslFunc.function = this->CCallback_Rateintegrandt;
    gslFunc.params = Callback_params;
    
    //gsl_integration_qags(&gslFunc, t_min, t_max, eps, 1e-4, maxInteration, gsl_workSpace, &gslresult, &gslerror);
    int gslQAGkey = 2;
    status = gsl_integration_qag(&gslFunc, t_min, t_max, eps, 1e-5, maxInteration, gslQAGkey, gsl_workSpace, &gslresult, &gslerror);


    gsl_integration_workspace_free(gsl_workSpace);
    delete Callback_params;
    delete [] paramsPtr;

    return(gslresult);
}

double HG_1to3_decay::Rateintegrandt(double t, void *params)
{
    double *par = (double*)params;
    double rateType = par[0];
    double Temp = par[1];
    double Eq = par[2];
    double s = par[3];

    double E1_min;
    double u = - s - t + m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
    E1_min = Eq*m[0]*m[0]/(m[0]*m[0] - u) + (m[0]*m[0] - u)/4/Eq;

    double *paramsPtr = new double [5];
    paramsPtr[0] = rateType;
    paramsPtr[1] = Temp;
    paramsPtr[2] = Eq;
    paramsPtr[3] = s;
    paramsPtr[4] = t;
    CCallbackHolder *Callback_params = new CCallbackHolder;
    Callback_params->clsPtr = this;
    Callback_params->params = paramsPtr;
    int maxInteration = 1000;
    gsl_integration_workspace *gsl_workSpace = gsl_integration_workspace_alloc(maxInteration);
    double gslresult, gslerror;
    int status;
    gsl_function gslFunc;
    gslFunc.function = this->CCallback_RateintegrandE1;
    gslFunc.params = Callback_params;
    
    status = gsl_integration_qagiu(&gslFunc, E1_min, eps, 1e-5, maxInteration, gsl_workSpace, &gslresult, &gslerror);

    gsl_integration_workspace_free(gsl_workSpace);
    delete Callback_params;
    delete [] paramsPtr;

    double matrixElementsSq = Matrix_elements_sq(s, t);

    return(gslresult*matrixElementsSq);
}

double HG_1to3_decay::RateintegrandE1(double E1, void *params)
{
    double result;
    double *par = (double*)params;
    double rateType = par[0];
    double Temp = par[1];
    double Eq = par[2];
    double s = par[3];
    double t = par[4];

    double E2_min;
    double E2_max;
    double min_1 = Eq*m[1]*m[1]/(t - m[1]*m[1]) + (t - m[1]*m[1])/4/Eq;

    double a = - (s + t - m[1]*m[1] - m[2]*m[2])*(s + t - m[1]*m[1] - m[2]*m[2]);
    double b = - Eq*((s + t - m[1]*m[1] - m[2]*m[2])*(s - m[0]*m[0] - m[1]*m[1]) 
               - 2*m[0]*m[0]*(m[1]*m[1] - t)) - E1*(m[1]*m[1] - t)
               *(s + t - m[1]*m[1] - m[2]*m[2]);
    double c = - (t - m[1]*m[1])*(t - m[1]*m[1])*E1*E1
               - 2*Eq*(2*m[1]*m[1]*(s + t - m[1]*m[1] - m[2]*m[2]) 
               - (m[1]*m[1] - t)*(s - m[0]*m[0] - m[1]*m[1]))*E1
               + 4*Eq*Eq*m[0]*m[0]*m[1]*m[1] + m[1]*m[1]*(s + t - m[1]*m[1] 
               - m[2]*m[2])*(s + t - m[1]*m[1] - m[2]*m[2]) + m[0]*m[0]
               *(m[1]*m[1] -t)*(m[1]*m[1] -t)
               - Eq*Eq*(s - m[0]*m[0] - m[1]*m[1])*(s - m[0]*m[0] - m[1]*m[1])
               + (s - m[0]*m[0] - m[1]*m[1])*(t - m[1]*m[1])*(s + t - m[1]*m[1] 
               - m[2]*m[2]);

    
    if((b*b - a*c) >= 0) 
    {
       int integrandForm = 1;
       double min_2 = (-b + sqrt(b*b - a*c))/a;
       if(min_1 < min_2)
       {
          E2_min = min_2;
          integrandForm = 2;
       }
       else
          E2_min = min_1;
       E2_max = (-b - sqrt(b*b - a*c))/a;

       if(E2_max < E2_min) return(0.0);

       double *paramsPtr = new double [9];
       paramsPtr[0] = rateType;
       paramsPtr[1] = Temp;
       paramsPtr[2] = Eq;
       paramsPtr[3] = s;
       paramsPtr[4] = t;
       paramsPtr[5] = E1;
       paramsPtr[6] = integrandForm;
       paramsPtr[7] = min_2;
       paramsPtr[8] = a;
       CCallbackHolder *Callback_params = new CCallbackHolder;
       Callback_params->clsPtr = this;
       Callback_params->params = paramsPtr;
       int maxInteration = 1000;
       gsl_integration_workspace *gsl_workSpace = gsl_integration_workspace_alloc(maxInteration);
       double gslresult, gslerror;
       int status;
       gsl_function gslFunc;
       gslFunc.function = this->CCallback_RateintegrandE2;
       gslFunc.params = Callback_params;
       gsl_integration_qaws_table* gslQAWSptr;
       if(integrandForm == 2)
       {
          gslQAWSptr = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0, 0);
          status = gsl_integration_qaws(&gslFunc, E2_min+eps, E2_max-eps, gslQAWSptr, eps, 1e-5, maxInteration, gsl_workSpace, &gslresult, &gslerror);
       }
       else
       {
          gslQAWSptr = gsl_integration_qaws_table_alloc(0.0, -0.5, 0, 0);
          status = gsl_integration_qaws(&gslFunc, E2_min+eps, E2_max-eps, gslQAWSptr, eps, 1e-5, maxInteration, gsl_workSpace, &gslresult, &gslerror);
       }
       
       //gsl_integration_qags(&gslFunc, E1_min+eps, E2_max-eps, 0, 1e-4, maxInteration, gsl_workSpace, &gslresult, &gslerror);

       gsl_integration_qaws_table_free(gslQAWSptr);
       gsl_integration_workspace_free(gsl_workSpace);
       delete Callback_params;
       delete [] paramsPtr;

       result = gslresult;
    }
    else  // no kinematic phase space
    {
       result = 0.0e0;
    }
    return(result);
}

double HG_1to3_decay::RateintegrandE2(double E2, void *params)
{
    double *par = (double*)params;
    int rateType = (int) par[0];
    double Temp = par[1];
    double Eq = par[2];
    double s = par[3];
    double t = par[4];
    double E1 = par[5];
    int integrandForm = (int) par[6];
    double min_2 = par[7];
    double a = par[8];

    double mu1 = mu[0];
    double mu2 = mu[1];
    double mu3 = mu[2];
    double f0_E1 = Bose_distribution(E1, Temp, mu1);
    double f0_E2 = Bose_distribution(E2, Temp, mu2);
    double f0_E3 = Bose_distribution(E1 - E2 - Eq, Temp, mu3);
    //double common_factor = f0_E1*(1 + f0_E2)*(1 + f0_E3)/(sqrt(a*E2*E2 + 2*b*E2 + c) + eps);
    double common_factor;
    if(integrandForm == 2)
       common_factor = f0_E1*(1 + f0_E2)*(1 + f0_E3)/sqrt(-a);
    else
       common_factor = f0_E1*(1 + f0_E2)*(1 + f0_E3)/sqrt((-a)*(E2 - min_2));

 
    double result;
    if(rateType == 0)
       result = common_factor;
    else if(rateType == 1)
       result = common_factor*viscous_integrand(s, t, E1, E2, Eq, Temp, f0_E1, f0_E2, f0_E3);
    else
    {
       if(bulk_deltaf_kind == 0)
       {
           double* bulkvis_B0 = new double [3];
           double* bulkvis_D0 = new double [3];
           double* bulkvis_E0 = new double [3];
           get_bulkvis_coefficients_14moment(Temp, bulkvis_B0, bulkvis_D0, bulkvis_E0);
           result = common_factor*bulkvis_integrand_14moment(E1, E2, Eq, f0_E1, f0_E2, f0_E3, bulkvis_B0, bulkvis_D0, bulkvis_E0);
           delete[] bulkvis_B0;
           delete[] bulkvis_D0;
           delete[] bulkvis_E0;
       }
       else if(bulk_deltaf_kind == 1)
       {
           double bulkvis_Cbulk, bulkvis_e2;
           get_bulkvis_coefficients_relaxation(Temp, &bulkvis_Cbulk, &bulkvis_e2);
           result = common_factor*bulkvis_integrand_relaxation(Temp, E1, E2, Eq, f0_E1, f0_E2, f0_E3, bulkvis_Cbulk, bulkvis_e2);
       }
    }

    return(result);
}

double HG_1to3_decay::Matrix_elements_sq(double s, double t)
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

double HG_1to3_decay::Matrix_elements_sq_C1(double s, double t)
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

double HG_1to3_decay::Matrix_elements_sq_C1_omega(double s, double t)
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

double HG_1to3_decay::Matrix_elements_sq_C2(double s, double t)
{
  //pi + pi -> rho + gamma
  double isospin_factor_C2p1 = 1.0;
  double isospin_factor_C2p2 = 2.0;
  double result = (C2p1(s,t)*isospin_factor_C2p1
                   + C2p2(s,t)*isospin_factor_C2p2);
  return (result);  
}

double HG_1to3_decay::Matrix_elements_sq_C3(double s, double t)
{
  //rho -> pi + pi + gamma
  double isospin_factor_C3p1 = 1.0;
  double isospin_factor_C3p2 = 2.0;
  double result = (C3p1(s,t)*isospin_factor_C3p1
                   + C3p2(s,t)*isospin_factor_C3p2);
  return (result);  
}

double HG_1to3_decay::Matrix_elements_sq_C4(double s, double t)
{
  //pi + Kstar -> K + gamma
  double isospin_factor_C4p1 = 4.0;
  double isospin_factor_C4p2 = 4.0;
  double result = (C4p1(s,t)*isospin_factor_C4p1
                   + C4p2(s,t)*isospin_factor_C4p2);
  return (result);
}

double HG_1to3_decay::Matrix_elements_sq_C5(double s, double t)
{
  //pi + K -> Kstar + gamma
  double isospin_factor_C5p1 = 4.0;
  double isospin_factor_C5p2 = 4.0;
  double result = (C5p1(s,t)*isospin_factor_C5p1
                   + C5p2(s,t)*isospin_factor_C5p2);
  return (result);
}

double HG_1to3_decay::Matrix_elements_sq_C6(double s, double t)
{
  //rho + K -> K + gamma
  double isospin_factor_C6p1 = 4.0;
  double isospin_factor_C6p2 = 4.0;
  double result = (C6p1(s,t)*isospin_factor_C6p1
                   + C6p2(s,t)*isospin_factor_C6p2);
  return (result);
}

double HG_1to3_decay::Matrix_elements_sq_C7(double s, double t)
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
double HG_1to3_decay::C1p1(double s, double t)
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

double HG_1to3_decay::C1p2(double s, double t)
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

double HG_1to3_decay::C1p3(double s, double t)
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

double HG_1to3_decay::C1p4(double s, double t)
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

double HG_1to3_decay::C1p5(double s, double t)
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

double HG_1to3_decay::C1p6(double s, double t)
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

double HG_1to3_decay::C2p1(double s, double t)
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

double HG_1to3_decay::C2p2(double s, double t)
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

double HG_1to3_decay::C3p1(double s, double t)
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

double HG_1to3_decay::C3p2(double s, double t)
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

double HG_1to3_decay::C4p1(double s, double t)
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

double HG_1to3_decay::C4p2(double s, double t)
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

double HG_1to3_decay::C5p1(double s, double t)
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

double HG_1to3_decay::C5p2(double s, double t)
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

double HG_1to3_decay::C6p1(double s, double t)
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

double HG_1to3_decay::C6p2(double s, double t)
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

double HG_1to3_decay::C7p1(double s, double t)
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

double HG_1to3_decay::C7p2(double s, double t)
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

double HG_1to3_decay::Power(double x, int a)
{
    return(pow(x,a));
}
