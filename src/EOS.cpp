#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include "Arsenal.h"
#include "Physicalconstants.h"
#include "EOS.h"
using namespace std;

EOS::EOS(int EOS_kind_in)
{
   EOS_kind = EOS_kind_in;
   if(EOS_kind == 0)
       EOS_Table_ptr = new Table("chemical_potential_tb/s95p/s95p-v1/EOS_PST.dat");
   else if(EOS_kind == 1)
       EOS_Table_ptr = new Table("chemical_potential_tb/s95p/s95p-PCE150-v1/EOS_PST.dat");
   else if(EOS_kind == 2)
       EOS_Table_ptr = new Table("chemical_potential_tb/s95p/s95p-PCE165-v0/EOS_PST.dat");
}

EOS::~EOS()
{
   delete EOS_Table_ptr;
}

double EOS::get_energy_density_from_temperature(double T)
{
   double ed = EOS_Table_ptr->interp(4, 1, T, 2);
   return(ed);
}

double EOS::get_pressure_from_temperature(double T)
{
   double pressure = EOS_Table_ptr->interp(4, 2, T, 2);
   return(pressure);
}

double EOS::get_cs2_from_temperature(double T)
{
   double ed = EOS_Table_ptr->interp(4, 1, T, 2);
   double de = ed/50.;
   double pLeft = EOS_Table_ptr->interp(1, 2, ed - de/2., 2);
   double pRight = EOS_Table_ptr->interp(1, 2, ed + de/2., 2);
   double cs2 = (pRight - pLeft)/de;
   return(cs2);
}
