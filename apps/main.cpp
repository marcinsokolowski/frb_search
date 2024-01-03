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
#include "dedisp_search.h"

string gOutputDir="candidates/";
string gImageFileNameTemplate="wsclean_%d_timeindex%03d-%04d-I-dirty.fits";

int gCoarseChannels=24;
int gTimeSteps=-1;    // 592;
int gStartTimeIndex=0; 
int gObsID=1192530256;
double gThresholdInSigma=10;
int gBorder=20;
mystring gSubDirTemplate="wsclean_timeindex%03d";

extern int gDebugLocal;

void usage()
{
   printf("frb_search -l MIN_DM -m MAX_DM -s DM_STEP -t TIME_STEPS -n THRESHOLD_IN_SIGMA -b BORDER -a ENABLE_ALGO_TYPE -d OUTDIR -o OBSID -f FITS_FILENAME_TEMPLATE -r MAX_ALLOWED_RMS -D\n");
   printf("\t-S START_TIME_INDEX : default = %d\n",gStartTimeIndex);
   printf("\t-t TIME_STEPS : set <=0 to analyse all available timesteps\n");
   printf("\t-a ENABLE_ALGO_TYPE : 1 - brute force, 2 - dynamic spectrum based also, 3 - both\n");
   printf("\t-o OBSID : obsid to be analysed [ default %d ]\n",gObsID);
   printf("\t-c CUT_THRESHOLD_IN_SIGMAS : threshold to cut values in dynamic spectrum in the dynaspec search [default %.2f]\n",CDedispSearch::m_ThresholdToCutInSigmas);   
   printf("\t-f FITS_FILENAME_TEMPLATE  : must contain %%d_timeindex%%03d-%%04d in the same order [default %s]\n",gImageFileNameTemplate.c_str());
   printf("\t-v : increase debug level by 1 [default debug level = %d]\n",gDebugLocal);
   printf("\t-r MAX_ALLOWED_RMS : default = %.4f\n",CDedispSearch::m_MaxRMSOnSingle);
   printf("\t-R EXCLUDE_REF_SOURCES : exclude reference sources [default %d]\n",CDedispSearch::m_bRejectRefSources);
   printf("\t-T Subdir template [default %s]\n",gSubDirTemplate.c_str());
   printf("\t-x Skip N first and last images [default %d]\n",CDedispSearch::m_SkipTimeSteps);
   printf("\t-D : enables de-dispersion debugging -> saves FITS files of pre-dedispersion images 0.5sec/1.28MHz [default %d]\n",CDedispSearch::m_bDebugDedispersion);
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ho:l:m:s:t:n:b:a:d:c:f:vr:R:S:T:x:D";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'a':
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
            break;

          case 'd':
            if( optarg ){
               gOutputDir = optarg;
            }
            break;

          case 'D':
            CDedispSearch::m_bDebugDedispersion = true;
            break;

          case 'f':
            if( optarg ){
               gImageFileNameTemplate = optarg;
            }
            break;

          case 'o':
            if( optarg ){
               gObsID = atol( optarg );
            }
            break;

         case 'n':
            if( optarg ){
               gThresholdInSigma = atof( optarg );
            }
            break;

         case 'c':
            if( optarg ){
               CDedispSearch::m_ThresholdToCutInSigmas = atof( optarg );
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

         case 'r':
            if( optarg ){
              CDedispSearch::m_MaxRMSOnSingle = atof( optarg );
            }
            break;

         case 'R':
            if( optarg ){
               CDedispSearch::m_bRejectRefSources = atol( optarg );
            }
            break;

         case 'S':
            if( optarg ){
               gStartTimeIndex = atol( optarg );
            }
            break;

         case 'T':
            if( optarg ){
               gSubDirTemplate = optarg;
            }
            break;

         case 'x':
            if( optarg ){
               CDedispSearch::m_SkipTimeSteps = atol( optarg );
            }
            break;

         case 'v':
            gDebugLocal++;
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

   CDedispSearch::m_ResultsPath = gOutputDir.c_str();
   
   // checks :
   if( fabs(CDedispSearch::m_StepDM ) < 0.00000000001 ){
      printf("ERROR : too small DM step - the loop might become infinite -> specify larger value in -s option!!!!\n");
      exit(-1);
      
   }
}
 

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("Start time index = %d\n",gStartTimeIndex);
   printf("Skip first images = %d\n",CDedispSearch::m_SkipTimeSteps);
   printf("ObsID     = %d\n",gObsID);   
   printf("DM range  = %.2f - %.2f\n",CDedispSearch::m_MinDM,CDedispSearch::m_MaxDM);
   printf("DM step   = %.2f\n",CDedispSearch::m_StepDM);
   printf("Max allowed RMS = %.4f\n",CDedispSearch::m_MaxRMSOnSingle);
   printf("Timesteps = %d\n",gTimeSteps);
   printf("Border    = %d\n",gBorder);
   printf("Algorithms :\n");
   printf("\tBrute force = %d\n",CDedispSearch::m_bRunBruteForceAlgo);
   printf("\t\tParameters:\n");
   printf("\tDynaapsec   = %d\n",CDedispSearch::m_bRunDynSpecAlgo);
   printf("\t\tParameters:\n");
   printf("\t\tCutInSigma = %.2f sigmas\n",CDedispSearch::m_ThresholdToCutInSigmas);
   printf("Outdir    = %s\n",gOutputDir.c_str());
   printf("Debug level = %d\n",gDebugLocal);
   printf("\t\tLog files used in de-dispersion = %d\n",CDedispSearch::m_bDebugDedispersion);
   printf("Exclude ref. sources = %d\n",CDedispSearch::m_bRejectRefSources);
   printf("Subdir template = %s\n",gSubDirTemplate.c_str());
   printf("#####################################\n");   
}



int main(int argc,char* argv[])
{
  time_t start_time = get_dttm();
  if( (argc>=2 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--h")==0)) ){
     usage();
  }

//  fits_flagged=argv[1];
//  fits_out = argv[2];

  parse_cmdline(argc+0,argv+0);
  print_parameters();

  CDedispSearch frb_search( gObsID, gCoarseChannels, gTimeSteps, gStartTimeIndex );
  int n_images = frb_search.m_MWADataCube.Read( gSubDirTemplate.c_str() , gImageFileNameTemplate.c_str() );
  int events_count = frb_search.Run( gThresholdInSigma, gBorder, 6 );
  
//  CMWADataCube data_cube( gObsID, gCoarseChannels, gTimeSteps );
//  int n_images = data_cube.Read();


   time_t end_time = get_dttm();
   printf("FRB finder finished at unixtime = %d, took %d seconds\n",(int)end_time,(int)(end_time-start_time));
}
