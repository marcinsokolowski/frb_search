#include <stdio.h>
#include <stdlib.h>
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
// #include "dedisp_search.h"

string gOutputDir="./";
string gInputFits="chan_204_20210224T20260800_XX.fits";

void usage()
{
   printf("read_data_cube FITS_FILE\n"); // -l MIN_DM -m MAX_DM -s DM_STEP -t TIME_STEPS -n THRESHOLD_IN_SIGMA -b BORDER -c Threshold_on_SNR\n");
//   printf("-c Threshold_on_SNR : threshold for SNR [default = %.2f]\n",CDedispSearch::m_SNRThreshold);
//   printf("-x Threshold_for_cleaning : threshold to ignore pixels with lower value [default %.2f]\n",CDedispSearch::m_SNRThreshold);
//   printf("-n THRESHOLD_IN_SIGMA : threshold to cut pixels in the dynamic spectrum [default %.2f]\n",gThresholdInSigma);
//   printf("-d DEBUG_LEVEL\n");
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ho:l:m:s:t:n:b:c:x:d:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'o':
            if( optarg ){
               gOutputDir = optarg;
            }
            break;


/*         case 'd':
            if( optarg ){
              gDebugLocal = atol( optarg );
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
   
//   string szFullPath=gOutputDir.c_str();
//   szFullPath += "/";
//   szFullPath += fits_out.c_str();
//   fits_out = szFullPath.c_str();
}
 

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("Input FITS = %s\n",gInputFits.c_str());
   printf("#####################################\n");   
}



int main(int argc,char* argv[])
{
  time_t start_time = get_dttm();
  if( (argc<2 || ( argc>=2 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--h")==0))) ){
     usage();
  }

  gInputFits = argv[1];

//  fits_flagged=argv[1];
//  fits_out = argv[2];

  parse_cmdline(argc+0,argv+0);
  print_parameters();


   const char* fits_file = argv[1]; 
   fitsfile *fp;
   int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
   int atype, btype, anaxis, bnaxis, check = 1, ii, op;
   long npixels = 1, firstpix[3] = {1,1,1}, ntodo;
//   long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1};
   double *apix, *bpix, value;
   int image2=1;
 
   printf("INFO : reading FITS cube %s\n",fits_file);
   fits_open_file( &fp, fits_file, READONLY, &status ); /* open input images */
   if (status) {
      fits_report_error(stderr, status); /* print error message */
      return(status);
   }

   fits_get_img_dim(fp, &anaxis, &status);  /* read dimensions */
   printf("INFO : image  %s has %d axis\n",fits_file,anaxis);
   
   long* anaxes = new long[anaxis];
   fits_get_img_size(fp, anaxis, anaxes, &status);
   if (status) {
      fits_report_error(stderr, status); /* print error message */
      return(status);
   }
   mystring szDim;
   for(int i=0;i<anaxis;i++){
      char szTmp[32];
      if ( i < (anaxis-1) ){
         sprintf(szTmp,"%d x ",anaxes[i]);
      }else{
         sprintf(szTmp,"%d",anaxes[i]);
      }
      szDim += szTmp;
   }
   printf("INFO : image dimensions = %s\n",szDim.c_str());

   printf("INFO : reading FITS header\n");
   CBgFits fits;   
   fits.ReadFits( gInputFits.c_str(), 0, 0 );
   printf("INFO : fits image %d x %d , freq_start = %.4f MHz\n",fits.GetXSize(),fits.GetYSize(),fits.start_freq);
   
   // test reading all channels :
   int sizeXY = fits.GetXSize()*fits.GetYSize();
   float* data = new float[sizeXY];
   long* firstpixel = new long[anaxis];
   int image_type = TFLOAT;
   firstpixel[0] = 1; // fits.GetXSize()
   firstpixel[1] = 1; // fits.GetYSize()
   firstpixel[3] = 1;
   
   
   int n_channels = anaxes[2];
   for(int channel=0;channel<n_channels;channel++){
      firstpixel[2] = channel + 1; // in FORTRAN 1-indexed
      
      fits_read_pix( fp, image_type, firstpixel, sizeXY, NULL, data, NULL, &status);      
   }
   
   
   CMWADataCube cube;
   if( cube.ReadFitsCube( gInputFits.c_str() ) ){
     printf("ERROR : while reading file %s\n",gInputFits.c_str());
   }



   // destructor
   if( anaxes ){
      delete [] anaxes;
   }
   if( data ){
      delete [] data;
   }
   if( firstpixel ){
      delete [] firstpixel;
   }
}
