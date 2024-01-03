#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>


#include <bg_globals.h>
#include <bg_fits.h>
#include <bg_array.h>
#include <bg_bedlam.h>

#include <myfile.h>
#include <mystring.h>

#include <vector>
using namespace std;

#include "mwa_fits.h"
#include "pixel.h"
#include "dedisp_search.h"

string gInFits="test.fits";

void usage()
{
   printf("calcstat FITS_FILE\n");
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ha:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
/*         case 'a':
            if( optarg ){            
               int param_val = atol( optarg );
               if( param_val == 1 ){
                   CDedispSearch::m_bRunBruteForceAlgo = 1;
               }
               if( param_val == 2 ){
                   CDedispSearch::m_bRunDynSpecAlgo = 1;
               }
               if( param_val == 3 ){
                   CDedispSearch::m_bRunBruteForceAlgo = 1;
                   CDedispSearch::m_bRunDynSpecAlgo = 1;
               }
               
            }
            break;*/

         case 'h':
            usage();
            break;
            
         default:  
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
   
} 

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("Input fits file = %s\n",gInFits.c_str());
   printf("#####################################\n");   
}



int main(int argc,char* argv[])
{
  time_t start_time = get_dttm();
  if( (argc>=2 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--h")==0)) ){
     usage();
  }


  gInFits = argv[1];
  parse_cmdline(argc-1,argv+1);
  print_parameters();


  CMWAFits in_fits( gInFits.c_str() );
  if( in_fits.ReadFits( gInFits.c_str() ) ){
     printf("ERROR : could read fits file %s\n", gInFits.c_str() );
     exit(-1);
  }
  
  double mean,rms, min_val, max_val;
  CPixel max_pixel, min_pixel;
  CBgArray max_along_x, max_along_y;
  vector<int> max_along_x_x, max_along_y_y;
  
  in_fits.GetStat( mean, rms, min_val, max_val, max_pixel, min_pixel, max_along_x, max_along_x_x, max_along_y, max_along_y_y );
  printf("###########################################################################\n");
  printf("mean    = %.4f\n",mean);
  printf("rms     = %.4f\n",rms);
  printf("min_val = %.4f at pixel (%d,%d) [debug value = %.4f]\n",min_val,min_pixel.time,min_pixel.cc,min_pixel.flux);
  printf("max_val = %.4f at pixel (%d,%d) [debug value = %.4f]\n",max_val,max_pixel.time,max_pixel.cc,max_pixel.flux);
  printf("###########################################################################\n");
  
//  max_along_x.SaveToFile("max_along_y.txt");
  MyOFile out_file_max_along_x( "max_along_y.txt" );
  for(int y=0;y<in_fits.GetYSize();y++){
     out_file_max_along_x.Printf("%d %.4f %d\n",(y+100),max_along_x[y],max_along_x_x[y]); // y+in_fits.inttime
  }
  
   
}
