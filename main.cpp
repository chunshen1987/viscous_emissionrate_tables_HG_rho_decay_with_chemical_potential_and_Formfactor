#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Arsenal.h"
#include "Matrix_elements_sq.h"
#include "Phasespace_integrals.h"
#include "parameters.h"
#include "Stopwatch.h"
#include "chemical_potential.h"

using namespace std;


int main()
{
   Stopwatch sw; 
   sw.tic();

   Chemical_potential* chempotential_ptr;
   Chemical_potential chemicalpotenital;
   chempotential_ptr = &chemicalpotenital;
   chempotential_ptr->readin_chempotential_table("chemical_potential_tb/s95p/s95p-PCE150-v1/s95p-v1-PCE150_chemvsT.dat");
   chempotential_ptr->Set_chemical_potential_s95pv0PCE();

   int channel = 0;
   double m[3];
   string filename;
   
   //C.3
   filename = "rho_to_pion_pion_gamma";
   channel = 8;
   m[0] = mrho;
   m[1] = mpion;
   m[2] = mpion;
   Calculate_emissionrates(m, chempotential_ptr, channel, filename);

   sw.toc();
   cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
