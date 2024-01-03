#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>


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
int gOBSID=1217495184;
double gThresholdInSigma=10.00;
double gThresholdInSigmaForSNR=10.00;
int gBorder=20;

string gInputFitsFile = "dyna_spec.fits";
string gOutFitsFile="series.fits";
string gOutCandidatesFile="cand.reg";
string gOutCountFitsFile = "count.fits";
string gOutAvgDynSpec = "dyna_spec_avg1.fits";

double gLogBelowDM = -1e6;
bool gIgnoreMaxValue = false;
bool gDoWeigthByRMS  = true;

bool gDivideByMapNew = true;

int gAverageNTimeSteps=1; // <=1 -> averaging in time turned off 

bool gZeros2NaNs = true;
double gZeroDistanceThreshold = 1e-9;

bool gAutoMetafits=false;
bool gTransposedFits=false; // true when dynamic spectra is "horizontal" time on X-axis and frequency on Y-axis

bool gAutoFlagging=false;
std::vector<string> gFlaggedRanges;
std::vector<CFlagRanges> gBadChannelRanges;

enum eDTSSearchType { eNoSearch=0, eNormalSearch=1 };
// do search or not :
eDTSSearchType gSearchType=eNormalSearch;

void usage()
{
   printf("find_sources_fits FITS_FILE_DTS\n");   
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "a:ho:l:m:s:t:n:N:b:d:o:L:Xw:rAD:TxF:fS:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
/*         case 'a':
            if( optarg ){
               gAverageNTimeSteps = atol( optarg );
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
   printf("Input fits file  = %s\n",gInputFitsFile.c_str());
   printf("Output fits file = %s\n",gOutFitsFile.c_str()); 
   printf("Output reg file  = %s\n",gOutCandidatesFile.c_str());
//   printf("OBSID     = %d\n",gOBSID);  
   printf("#####################################\n");   
}



int main(int argc,char* argv[])
{
  // initialising defaults for this particular program :
//  CDedispSearch::m_MinDM = 100.00;
//  CDedispSearch::m_MaxDM = 2000.00;

  time_t start_time = get_dttm();
  if( (argc<2 || ( argc>=2 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--h")==0))) ){
     usage();
  }

  gInputFitsFile = argv[1];
  
  parse_cmdline(argc-1,argv+1);
  print_parameters();


  // const char* change_ext(const char* name,const char* new_ext,string& out);  
  string gOutRmsMapFile, gOutSnrSeriesFile, gOutSnrSeriesRegFile, gOutRmsMapNewFile;
  // change_ext( gOutFitsFile.c_str(), "reg", gOutCandidatesFile );
  change_ext( gOutFitsFile.c_str(), "reg", gOutSnrSeriesRegFile ); // true ?

  CMWAFits dedispersed_series( gInputFitsFile.c_str() );
  dedispersed_series.ReadFits( gInputFitsFile.c_str(), 0, 1, 0, gTransposedFits );
  
  //    int FindSourcesSNR( double threshold_snr=5.00, int border = 20, const char* szOutRegFile="sources_snr.reg", int minX=0, int maxX=-1, CBgFits* pCountMap=NULL , bool bSavePhysical=false);
  dedispersed_series.FindSourcesSNR( gThresholdInSigmaForSNR, 0, gOutSnrSeriesRegFile.c_str(), 0, -1, NULL, true ); // ,  start_time_index+2, end_time_index-2, &count_map );
  
}
