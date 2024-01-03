#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <string>

using namespace std;

int main()
{
/*   string file1 = "chan_204_20210224T20261950_YY.fits";
   char file[1024];
   strcpy(file,file1.c_str());
   
   char* ptr = (char*)strstr(file,"YY");
   if( ptr ){
      ptr[0] = '_';
      ptr[1] = '_';
      
      printf("Output = |%s|\n",file);
   }*/
   
   double val = 1.00 + 0.00;
   
   printf("isnan(%.8f) = %d\n",val,isnan(val));
   
   if( isnan(val) == 0 ){
      printf("Value = %.8f is NOT NaN\n",val);
   }else{
      printf("Value = %.8f is NaN\n",val);
   }

   if( isinf(val) == 0 ){
      printf("Value = %.8f is NOT INF\n",val);
   }else{
      printf("Value = %.8f is INF\n",val);
   }
}