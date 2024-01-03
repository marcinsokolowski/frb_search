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

string gInputFitsFile = "dyna_spec.fits";
string gOutFitsFile="series.fits";
double gDM=100.00;
double gValue=5.00;
int    gOBSID=1217495184;
int    gStartTimeIndex=200;

void usage()
{
   printf("inject_frb FITS_FILE_DYNA-SPEC OUTPUT_FITSFILE -d DM -v VALUE -o OBSID -s START_TIMEINDEX\n");
   printf("-s START_TIMEINDEX : start time index [defualt %d]\n",gStartTimeIndex);   
//   printf("-v DEBUG_LEVEL\n");
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "hd:v:o:s:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'd':
            if( optarg ){
              gDM = atof( optarg );
            }
            break;

         case 'v':
            if( optarg ){
              gValue = atof( optarg );
            }
            break;

         case 'o':
            if( optarg ){
               gOBSID = atol( optarg );
            }
            break;

         case 's':
            if( optarg ){
               gStartTimeIndex = atol( optarg );
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

   CDedispSearch::m_MinDM = gDM;
   CDedispSearch::m_MaxDM = gDM;
   
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
   printf("Input dynamic spectrum file  = %s\n",gInputFitsFile.c_str());
   printf("OBSID = %d\n",gOBSID);
   printf("DM    = %.2f\n",gDM);
   printf("VALUE = %.2f\n",gValue);   
   printf("START_TIMEINDEX = %d\n",gStartTimeIndex);
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
  
  if( argc>=3 ){
     gOutFitsFile = argv[2];
  }

//  fits_flagged=argv[1];
//  fits_out = argv[2];

  parse_cmdline(argc-2,argv+2);
  print_parameters();
  

  CMWAFits dyna_spec( gInputFitsFile.c_str() );
  dyna_spec.ReadFits( gInputFitsFile.c_str() );

  CDedispSearch dedisp( gOBSID, dyna_spec.GetYSize(), dyna_spec.GetXSize() );
//  dedisp.m_MWADataCube.m_StartUnixTime = 
  dedisp.m_MWADataCube.m_Timesteps = dyna_spec.GetXSize();
  dedisp.m_MWADataCube.m_Channels  = dyna_spec.GetYSize();

  dedisp.m_MWADataCube.ReadMetaData( gOBSID );

  double delta_freq = 0.00;
  std::vector<double> freq_list;
  if( dedisp.m_MWADataCube.get_freq_list( freq_list ) != dedisp.m_MWADataCube.m_CoarseChannels.size() ||  dedisp.m_MWADataCube.m_CoarseChannels.size() <= 0 ){
     printf("ERROR : could not create the list of observing frequencies\n");
     exit(-1);
  }


  double inttime = (dedisp.m_MWADataCube.m_Metafits)->inttime;
  printf("inttime = %.2 sec\n",inttime);  
  dedisp.CalcTemplatesExact( dyna_spec );

  // insert FRB :
  dedisp.InjectFRB( dyna_spec, gStartTimeIndex, gValue );
    
  dyna_spec.WriteFits( gOutFitsFile.c_str() );

}
