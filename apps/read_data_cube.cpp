#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include <bg_globals.h>
#include <bg_fits.h>
#include <bg_array.h>
#include <bg_bedlam.h>
#include <libnova_interface.h>

#include <myfile.h>
#include <mystring.h>

#include <vector>
using namespace std;

#include "mwa_fits.h"
// #include "dedisp_search.h"

string gOutputDir="dynamic_spectra/";
string gInputFits="chan_204_20210224T20260800_XX.fits";
int gFreqChannel = 204;
double gFreqChannelDist = (400.00/512.00);
double gBW_MHz = (400.00/512.00)*(32.00/27.00);
double gInttime = 0.100; // in seconds 100 or 10 ms 
int gX = 89;
int gY = 89;
int gRadius=-1;

string gListYY;

bool bSaveRescaledFits=false;
int  bSaveImagesOnly=0;

void usage()
{
   printf("read_data_cube FITS_FILE FREQ_CHANNEL\n"); // -l MIN_DM -m MAX_DM -s DM_STEP -t TIME_STEPS -n THRESHOLD_IN_SIGMA -b BORDER -c Threshold_on_SNR\n");
   printf("-x X : position X of pixel to dump [default %d]\n",gX);
   printf("-y Y : position Y of pixel to dump [default %d]\n",gY);
   printf("-r RADIUS : radius in pixels around the pixels specified by -x and -y (default center) to dump dynamic spectra and filerbank files [default <0 -> just dump a single pixel]\n");
   printf("-i INTTIME : integration time [default %.3f seconds]\n",gInttime);
   printf("-Y LIST_FITS_Y : list of YY images to calculate Stokes I = (XX+YY)/2\n");
   printf("-R : save rescaled fits files [default %d]\n",bSaveRescaledFits);
   printf("-S : save images only [default =0 which is disabled], set to >0 to save specified number of timesteps or -1 to save all timesteps in data cube\n");
//   printf("-c Threshold_on_SNR : threshold for SNR [default = %.2f]\n",CDedispSearch::m_SNRThreshold);
//   printf("-x Threshold_for_cleaning : threshold to ignore pixels with lower value [default %.2f]\n",CDedispSearch::m_SNRThreshold);
//   printf("-n THRESHOLD_IN_SIGMA : threshold to cut pixels in the dynamic spectrum [default %.2f]\n",gThresholdInSigma);
//   printf("-d DEBUG_LEVEL\n");
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ho:x:y:r:i:Y:RS:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'i':
            if( optarg ){
               gInttime = atof( optarg );
            }
            break;

         case 'o':
            if( optarg ){
               gOutputDir = optarg;
            }
            break;

         case 'x':
            if( optarg ){
               gX = atol( optarg );
            }
            break;

         case 'y':
            if( optarg ){
               gY = atol( optarg );
            }
            break;

         case 'r':
            if( optarg ){
               gRadius = atol( optarg );
            }
            break;

         case 'R':
            bSaveRescaledFits = true;
            break;

         case 'Y':
            if( optarg ){
               gListYY = optarg;
            }
            break;

         case 'S':
            bSaveImagesOnly = atol( optarg );
            break;


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

   if( strlen(gOutputDir.c_str()) && strcmp(gOutputDir.c_str(),"./") ){
      MyFile::CreateDir( gOutputDir.c_str() );
   }
}
 

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("Input list = %s\n",gInputFits.c_str());
   printf("List of YY images = %s\n",gListYY.c_str());
   printf("Dump dynamic spectrum of a pixel at (x,y) = (%d,%d)\n",gX,gY);
   printf("Radius = %d\n",gRadius);
   printf("Integration time = %.3f [sec]\n",gInttime);
   printf("Re-scale FITS files = %d\n",bSaveRescaledFits);
   printf("Save images only    = %d\n",bSaveImagesOnly);
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
  string filfile_out = argv[2];
//  string spec_file_out = "spectrum.txt";
  
  parse_cmdline(argc-2,argv+2);
  print_parameters();


  const char* fits_file = argv[1]; 

  double freq_center = gFreqChannel*gFreqChannelDist;
  double freq_lower_end = freq_center - gBW_MHz/2.00;
  double freq_upper_end = freq_center + gBW_MHz/2.00;

  CMWADataCube cube;   
  if( cube.ReadFitsCube( gInputFits.c_str(), freq_center, freq_lower_end, freq_upper_end ) ){
     printf("ERROR : while reading data cube from FITS files in the list  %s\n",gInputFits.c_str());
     exit(-1);
  }
  
  if( bSaveImagesOnly != 0 ){
     cube.SaveChannelFitsFiles( bSaveImagesOnly );
     exit(0);
  }
  
  if( strlen(gListYY.c_str()) > 0 ){
     CMWADataCube cube_yy;
     
     if( cube_yy.ReadFitsCube( gListYY.c_str(), freq_center, freq_lower_end, freq_upper_end ) ){
         printf("ERROR : while reading data cube from FITS files in the list  %s\n",gListYY.c_str());
         exit(-1);
     }
     
     printf("INFO : calculating Stokes I = (XX+YY)/2 - using images in list files %s and %s\n",gInputFits.c_str(),gListYY.c_str());
     cube.Avg( cube_yy );
  }
  
  double foff = -(freq_upper_end - freq_lower_end) / cube.m_Channels;
  double tstart = ux2mjd( int(cube.m_StartUnixTime) , (cube.m_StartUnixTime-int(cube.m_StartUnixTime))*1000000.00 );
  printf("DEBUG : freq_upper_end = %.4f MHz, foff = %.4f MHz, tstart = %.4f MJD\n",freq_upper_end, foff, tstart);
  
  if( gRadius > 0 ){
     // dump pixels in a radius around (x,y) :
     printf("INFO : generating dynamic spectra and filterbank files in radius %d pixels around (%d,%d)\n",gRadius,gX,gY);
     
     for(int yy=(gY-gRadius-1);yy<=(gY+gRadius+1);yy++){
        for(int xx=(gX-gRadius-1);xx<=(gX+gRadius+1);xx++){
           if( yy>=0 && yy<cube.GetYSize() && xx>=0 && xx<cube.GetXSize() ){
              double dist = sqrt( (xx-gX)*(xx-gX) + (yy-gY)*(yy-gY) );
              if( dist <= gRadius ){
                 printf("PROGRESS : dumping pixel (%d,%d)\n",xx,yy);
              
                 char szOutFilFile[1024];
                 sprintf(szOutFilFile,"%s/dynaspec_%05d_%05d.fil",gOutputDir.c_str(),xx,yy);
                 filfile_out = szOutFilFile;
              
                 cube.generate_filfile( xx , yy, filfile_out.c_str(), freq_upper_end, foff, cube.m_StartUnixTime, gBW_MHz, gInttime, bSaveRescaledFits );
              }
           }
        }
     }
  }else{
     // just dump one pixel 
     cube.generate_filfile( gX , gY, filfile_out.c_str(), freq_upper_end, foff, tstart, gBW_MHz, gInttime, bSaveRescaledFits );
  }

}
