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
#include "dedisp_search.h"

string gOutputDir="./";

int gCoarseChannels=24;
int gTimeSteps=592;
int gObsID=1192530256;
double gThresholdInSigma=-100;
int gBorder=20;

void usage()
{
   printf("dynspec_test FITS_FILE_DYNA-SPEC -l MIN_DM -m MAX_DM -s DM_STEP -t TIME_STEPS -n THRESHOLD_IN_SIGMA -b BORDER -c Threshold_on_SNR\n");
   printf("-c Threshold_on_SNR : threshold for SNR [default = %.2f]\n",CDedispSearch::m_SNRThreshold);
   printf("-x Threshold_for_cleaning : threshold to ignore pixels with lower value [default %.2f]\n",CDedispSearch::m_SNRThreshold);
   printf("-n THRESHOLD_IN_SIGMA : threshold to cut pixels in the dynamic spectrum [default %.2f]\n",gThresholdInSigma);
   printf("-d DEBUG_LEVEL\n");
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

         case 'c':
            if( optarg ){
               CDedispSearch::m_SNRThreshold  = atof( optarg );
            }
            break;

         case 'n':
            if( optarg ){
               gThresholdInSigma = atof( optarg );
            }
            break;

         case 'l':
            if( optarg ){
               CDedispSearch::m_MinDM = atof( optarg );
            }
            break;

         case 'm':
            if( optarg ){
              CDedispSearch::m_MaxDM = atof( optarg );
            }
            break;

         case 's':
            if( optarg ){
              CDedispSearch::m_StepDM = atof( optarg );
            }
            break;

         case 't':
            if( optarg ){
              gTimeSteps = atol( optarg );
            }
            break;

         case 'b':
            if( optarg ){
              gBorder = atol( optarg );
            }
            break;

         case 'd':
            if( optarg ){
              gDebugLocal = atol( optarg );
            }
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
}
 

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("DM range  = %.2f - %.2f\n",CDedispSearch::m_MinDM,CDedispSearch::m_MaxDM);
   printf("DM step   = %.2f\n",CDedispSearch::m_StepDM);
   printf("Timesteps = %d\n",gTimeSteps);
   printf("Border    = %d\n",gBorder);
   printf("Threshold on SNR   = %.2f\n",CDedispSearch::m_SNRThreshold);
   printf("Threshold to clean = %.2f\n",gThresholdInSigma);
   printf("#####################################\n");   
}



int main(int argc,char* argv[])
{
  time_t start_time = get_dttm();
  if( (argc<2 || ( argc>=2 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--h")==0))) ){
     usage();
  }

  string fits_file = argv[1];

//  fits_flagged=argv[1];
//  fits_out = argv[2];

  parse_cmdline(argc+0,argv+0);
  print_parameters();


  CMWAFits dyna_spec( fits_file.c_str() );
  dyna_spec.ReadFits( fits_file.c_str() );

  CDedispSearch dedisp( 0, dyna_spec.GetYSize(), dyna_spec.GetXSize() );
//  dedisp.m_MWADataCube.m_StartUnixTime = 
  dedisp.m_MWADataCube.m_Timesteps = dyna_spec.GetXSize();
  dedisp.m_MWADataCube.m_Channels  = dyna_spec.GetYSize();
  double inttime = dyna_spec.GetIntTime();
  HeaderRecord* pFreqStart = dyna_spec.GetKeyword( "CRVAL2" );
  HeaderRecord* pDeltaFreq = dyna_spec.GetKeyword( "CDELT2" );
  double delta_freq = 0.00;
  std::vector<double> freq_list;
  
  if( pFreqStart && pDeltaFreq ){
      dyna_spec.start_freq = atof( pFreqStart->Value.c_str() );
      delta_freq =  atof( pDeltaFreq->Value.c_str() );
      dyna_spec.stop_freq  = dyna_spec.start_freq + dedisp.m_MWADataCube.m_Channels * delta_freq;      
      MWA_COARSE_CHANNEL = delta_freq;
      
      dedisp.m_MWADataCube.m_CoarseChannels.clear();
      for(int ch=0;ch<dedisp.m_MWADataCube.m_Channels;ch++){
         double freq_mhz = dyna_spec.start_freq + ch * delta_freq;
      
         dedisp.m_MWADataCube.m_CoarseChannels.push_back( ch );
         freq_list.push_back( freq_mhz );
      }
      
  }

      // test ch 121 :
/*      dedisp.m_MWADataCube.m_CoarseChannels.clear();
      for(int ch=0;ch<24;ch++){
         int cc = 109 + ch;
         double freq_mhz = cc*1.28;
         
         dedisp.m_MWADataCube.m_CoarseChannels.push_back( cc );
         freq_list.push_back( freq_mhz );
      }
      dyna_spec.start_freq = 109*1.28 - 0.64;
      delta_freq = 1.28;
      dyna_spec.stop_freq = 132*1.28 + 0.64;*/
      // end test 

  
  printf("Using inttime = %.6f [seconds]\n",inttime);
  printf("Freq range    = %.4f - %.4f [MHz]\n",dyna_spec.start_freq,dyna_spec.stop_freq);
  
//  CDedispSearch::m_MinDM = 10;
//  CDedispSearch::m_MaxDM = 150;

  dedisp.CalcTemplatesExact( inttime , dyna_spec.stop_freq , &freq_list );

  MyOFile candidates_file( "candidates.txt", "w" );
  candidates_file.Printf( LOG_FILE_HEADER );
  
  
  int min_good_pixels = 12; // 20;
  dedisp.AnalyseDynSpec( dyna_spec, 0, 0, gThresholdInSigma, &candidates_file, 1, min_good_pixels, &freq_list );

}
