#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include "Arsenal.h"
#include "parameters.h"
#include "chemical_potential.h"
using namespace std;

Chemical_potential::Chemical_potential()
{
   return;
}

Chemical_potential::~Chemical_potential()
{
   delete[] T;
   delete[] mu_pion;
   delete[] mu_K;
   delete EOS_Mu_Table_ptr;
   return;
}

void Chemical_potential::readin_chempotential_table(string filename)
{
    ostringstream filename_stream;
    filename_stream << filename;
    EOS_Mu_Table_ptr = new Table2D(filename_stream.str().c_str());

    return;
}

void Chemical_potential::Set_chemical_potential_s95pv0PCE()
{
    int Tb_sizeX = EOS_Mu_Table_ptr->getTbsizeX();
    Tb_length = Tb_sizeX;

    T = new double [Tb_sizeX];
    mu_pion = new double [Tb_sizeX];
    mu_K = new double [Tb_sizeX];

    for(int i=0; i<Tb_sizeX ; i++)
    {
       T[i] = EOS_Mu_Table_ptr->getTbdata(i, 0);
       mu_pion[i] = EOS_Mu_Table_ptr->getTbdata(i, 1);
       mu_K[i] = EOS_Mu_Table_ptr->getTbdata(i, 4);
    }

    return;
}

void Chemical_potential::Calculate_mu(double* Temperature, double* mu1, double* mu2, double* mu3, int npoint, int channel)
{
    switch(channel)
    {
       case 1:  //C1: pi + rho -> pi + gamma
       case 2:  //C1_omega: pi + rho -> omega -> pi + gamma
          interpolation1D_linear(T, mu_pion, Temperature, mu2, Tb_length, npoint);
          for(int i=0; i<npoint; i++)
          {
             mu1[i] = 2*mu2[i];
             mu3[i] = mu2[i];
          }
          break;
       case 3: //C2: pi + pi -> rho + gamma
          interpolation1D_linear(T, mu_pion, Temperature, mu1, Tb_length, npoint);
          for(int i=0; i<npoint; i++)
          {
             mu2[i] = mu1[i];
             mu3[i] = 2*mu1[i];
          }
          break;
       case 4: //C4: pi + Kstar -> K + gamma
          interpolation1D_linear(T, mu_pion, Temperature, mu2, Tb_length, npoint);
          interpolation1D_linear(T, mu_K, Temperature, mu3, Tb_length, npoint);
          for(int i=0; i<npoint; i++)
          {
             mu1[i] = mu2[i] + mu3[i];
          }
          break;
       case 5: //C5: pi + K -> Kstar + gamma
          interpolation1D_linear(T, mu_K, Temperature, mu1, Tb_length, npoint);
          interpolation1D_linear(T, mu_pion, Temperature, mu2, Tb_length, npoint);
          for(int i=0; i<npoint; i++)
          {
             mu3[i] = mu1[i] + mu2[i];
          }
          break;
       case 6: //C6: rho + K -> K + gamma
          interpolation1D_linear(T, mu_pion, Temperature, mu1, Tb_length, npoint);
          interpolation1D_linear(T, mu_K, Temperature, mu2, Tb_length, npoint);
          for(int i=0; i<npoint; i++)
          {
             mu1[i] = mu1[i]*2.;
             mu3[i] = mu2[i];
          }
          break;
       case 7: //C7: K + Kstar -> pi + gamma
          interpolation1D_linear(T, mu_K, Temperature, mu2, Tb_length, npoint);
          interpolation1D_linear(T, mu_pion, Temperature, mu3, Tb_length, npoint);
          for(int i=0; i<npoint; i++)
          {
             mu1[i] = mu2[i] + mu3[i];
          }
          break;
       case 8: //C3: rho -> pi + pi + gamma
          interpolation1D_linear(T, mu_pion, Temperature, mu2, Tb_length, npoint);
          for(int i=0; i<npoint; i++)
          {
             mu1[i] = 2.*mu2[i];
             mu3[i] = mu2[i];
          }
          break;
       default:
          cout << "calcualte chemical potential ERROR: can not find the corresponding channel! (channel =  " << channel << ")" << endl;
          exit(1);
    }
    return;
}

