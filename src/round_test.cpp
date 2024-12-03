#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double my_dm_index_round( double dm_index )
{
   int dm_index_int = int(dm_index);
   double frac = dm_index - dm_index_int;
   printf("DEBUG : frac = %.20f\n",frac);
   
   // <=0.5 -> int()
   if (frac <= 0.50000000001 ){
      printf("??? BELOW 0.5 -> YES\n");
      return dm_index_int;
   }
      
   return round(dm_index);
}


int main( int argc, char* argv[] )
{
   double m_MinDM = 0.00;
   double m_StepDM = 0.1;
   double dm = 0.1;
   if( argc >= 2 ){
      dm = atof( argv[1] );
   }
   
//   int out_dm_index = round( (dm-m_MinDM-m_StepDM/2.00) / m_StepDM);
   double value = (dm-m_MinDM-m_StepDM/2.00) / m_StepDM;
   printf("Rounding value = %.8f\n",value);
   int out_dm_index = my_dm_index_round( value );
   
   printf("out_dm_index = %d\n",out_dm_index);
}
