#ifndef CHEMICAL_POTENTIAL_H
#define CHEMICAL_POTENTIAL_H

#include "Table2D.h"

class Chemical_potential
{
   private: 
      int EOS_PCE_kind;
      int Tb_length;
      Table2D* EOS_Mu_Table_ptr;
      double* T;
      double* mu_pion;
      double* mu_K;
      
   public:
      Chemical_potential(int EOS_PCE_kind_in);
      ~Chemical_potential();

      void readin_chempotential_table(string filename);
      void Set_chemical_potential();

      double get_Tblength() {return(Tb_length);};
      double get_T(int i) {return(T[i]);};
      double get_mu_pion(int i) {return(mu_pion[i]);};
      double get_mu_K(int i) {return(mu_K[i]);};

      void Calculate_mu(double* Temperature, double* mu1, double* mu2, double* mu3, int npoint, int channel); 

};


#endif
