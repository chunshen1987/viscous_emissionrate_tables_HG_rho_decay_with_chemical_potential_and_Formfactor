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
#include "EOS.h"

using namespace std;


int main(int argc, char** argv)
{
   Stopwatch sw; 
   sw.tic();
   
   ParameterReader* paraRdr = new ParameterReader();
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);

   int EOS_PCE_kind = paraRdr->getVal("EOS_PCE_kind");

   Chemical_potential* chempotential_ptr = new Chemical_potential(EOS_PCE_kind);
   EOS* EOS_ptr = new EOS(EOS_PCE_kind);

   HG_1to3_decay HG1to3Rates(paraRdr, EOS_ptr);
   
   int channel = paraRdr->getVal("channel");
   HG1to3Rates.Calculate_emissionrates(chempotential_ptr, channel);

   sw.toc();
   cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
