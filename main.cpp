#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Arsenal.h"
#include "Physicalconstants.h"
#include "Stopwatch.h"
#include "chemical_potential.h"
#include "HG_1to3_decay.h"

using namespace std;


int main(int argc, char** argv)
{
   Stopwatch sw; 
   sw.tic();

   Chemical_potential* chempotential_ptr;
   Chemical_potential chemicalpotenital;
   chempotential_ptr = &chemicalpotenital;
   chempotential_ptr->readin_chempotential_table("chemical_potential_tb/s95p/s95p-PCE165-v0/s95p-v0-PCE165_chemvsT.dat");
   chempotential_ptr->Set_chemical_potential_s95pv0PCE();

   int channel = 0;
   double m[3];
   string filename;
   
   ParameterReader* paraRdr = new ParameterReader();
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);

   HG_1to3_decay HG1to3Rates(paraRdr);

   //C.3
   filename = "rho_to_pion_pion_gamma";
   channel = 8;
   HG1to3Rates.Calculate_emissionrates(chempotential_ptr, channel, filename);

   sw.toc();
   cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
