#ifndef EOS_H
#define EOS_H

#include "Table.h"

class EOS
{
   private: 
      int EOS_kind;
      Table* EOS_Table_ptr;
      
   public:
      EOS(int EOS_kind_in);
      ~EOS();
      double get_energy_density_from_temperature(double T);
      double get_pressure_from_temperature(double T);
      double get_cs2_from_temperature(double T);
};


#endif
