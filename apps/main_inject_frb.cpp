#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <bg_globals.h>
#include <bg_fits.h>
#include <bg_array.h>
#include <bg_bedlam.h>

#include <myfile.h>
#include <mystring.h>
#include <random.h>

#include <vector>
using namespace std;

#include "mwa_fits.h"
#include "dedisp_search.h"


string gOutputDir="./";

string gInputFitsFile = "dyna_spec.fits";
string gOutFitsFile="series.fits";
double gDM=100.00;
double gValue=-1.00;
double gSNR=-1;
int    gOBSID=1217495184;
int    gStartTimeIndex=-1; // <0 means random time

double gFrbFluence = 200.00;
double gIntTimeMiliseconds = 500.0; // 0.5 seconds 

vector<double> gDM_list;
vector<double> gFluenceList;

int gNumberOfFRBs = 1;

bool gAutoMetafits=false;

void usage()
{
   printf("inject_frb FITS_FILE_DYNA-SPEC OUTPUT_FITSFILE -d DM -v VALUE -o OBSID -s START_TIMEINDEX\n");
   printf("-s START_TIMEINDEX : start time index [defualt %d]\n",gStartTimeIndex);   
   printf("-f FLUENCE in [Jy ms] [default %.2f Jy ms], inserted value of flux density in Jy is calculated as FLUENCE/INTEGRATION_TIME_in_ms\n",gFrbFluence);
   printf("-i integration time in miliseconds [default %.1f ms = %.3f sec]\n",gIntTimeMiliseconds,(gIntTimeMiliseconds/1000.00));
   printf("-v VALUE : flux density of the inserted FRB pulse, overwrites option -f [default not set]\n");
   printf("-S SNR_VALUE : overwrites -v option and specifies what SNR should have the injected FRB signal\n");
   printf("-n NUMBER : number of FRBs to generated default = %d, if >1 -> WARNING: many things are still hardcoded and name of output FITS are automatically generated\n",gNumberOfFRBs);
   printf("-A : auto metafits fits\n");
   printf("-V : increase verbosity level\n");
//   printf("-v DEBUG_LEVEL\n");
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "hd:v:o:s:f:i:n:AVS:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'A':
            gAutoMetafits = true;
            break;
      
         case 'd':
            if( optarg ){
              gDM = atof( optarg );
            }
            break;

         case 'i':
            if( optarg ){
              gIntTimeMiliseconds = atof( optarg );
            }
            break;

         case 'f':
            if( optarg ){
              gFrbFluence = atof( optarg );
            }
            break;

         case 'v':
            if( optarg ){
              gValue = atof( optarg );
            }
            break;

         case 'S':
            if( optarg ){
              gSNR = atof( optarg );
            }
            break;

         case 'n':
            if( optarg ){
               gNumberOfFRBs = atol( optarg );
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

         case 'V':
            gBGPrintfLevel++;
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
   
   if( gValue < 0 ){
      gValue = gFrbFluence / gIntTimeMiliseconds;
   }
   
   if( gNumberOfFRBs > 1 ){
// 20191109 list :
/*      gDM_list.push_back(300);
      gDM_list.push_back(500);
      gDM_list.push_back(1000.00);
      gDM_list.push_back(2000.00);
      gDM_list.push_back(2500.00);
      gDM_list.push_back(2900.00);
      
      gFluenceList.push_back(  50.00 );
      gFluenceList.push_back( 100.00 );
      gFluenceList.push_back( 200.00 );
      gFluenceList.push_back( 300.00 );
      gFluenceList.push_back( 400.00 );
      gFluenceList.push_back( 500.00 ); // 1 Jy 
      gFluenceList.push_back( 1000.00 ); // 2 Jy 
      gFluenceList.push_back( 2000.00 ); // 4 Jy 
      gFluenceList.push_back( 3000.00 ); // 6 Jy 
      gFluenceList.push_back( 4000.00 ); // 8 Jy 
      gFluenceList.push_back( 5000.00 ); // 10 Jy 
*/
      
// systematic :
//      gFluenceList.push_back(  50.00 );
//      gFluenceList.push_back( 100.00 );
//      gFluenceList.push_back( 200.00 );
//      gFluenceList.push_back( 500.00 );
//      gFluenceList.push_back( 1000.00 );
//      gFluenceList.push_back( 2500.00 );      
      gFluenceList.push_back( gFrbFluence );     
      
      if( gDM_list.size() <= 0 ){
         gDM_list.push_back( gDM );
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
   printf("Input dynamic spectrum file  = %s\n",gInputFitsFile.c_str());
   printf("OBSID = %d\n",gOBSID);
   printf("DM    = %.2f\n",gDM);
   printf("Injected signal:\n");
   printf("\tVALUE = %.4f [Jy]\n",gValue);   
   printf("\tSNR   = %.4f\n",gSNR);
   printf("Fluence = %.2f [Jy ms]\n",gFrbFluence);
   printf("Integration time = %.2f [miliseconds] = %.4f [seconds]\n",gIntTimeMiliseconds,(gIntTimeMiliseconds/1000.00));
   printf("START_TIMEINDEX = %d\n",gStartTimeIndex);
   printf("Number of injections = %d\n",gNumberOfFRBs);
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
  
  change_ext( gInputFitsFile.c_str(), "_genfrb.fits", gOutFitsFile, true );
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

  if( gAutoMetafits ){
     printf("DEBUG : reading of metafits file is not required\n");
     dedisp.m_MWADataCube.GenerateMetaData( dyna_spec, dedisp, gOBSID );
  }else{
     dedisp.m_MWADataCube.ReadMetaData( gOBSID );
  }

  double delta_freq = 0.00;
  std::vector<double> freq_list;
  if( dedisp.m_MWADataCube.get_freq_list( freq_list ) != dedisp.m_MWADataCube.m_CoarseChannels.size() ||  dedisp.m_MWADataCube.m_CoarseChannels.size() <= 0 ){
     printf("ERROR : could not create the list of observing frequencies\n");
     exit(-1);
  }


  double inttime = (dedisp.m_MWADataCube.m_Metafits)->inttime;
  printf("inttime = %.2f sec\n",inttime);      
  
  if( gNumberOfFRBs <= 1 ){
     if( gNumberOfFRBs == 1 ){
        printf("Injecting just 1 FRB\n");
        dedisp.CalcTemplatesExact( dyna_spec );
        
        int start_index = gStartTimeIndex;
        if( start_index < 0 ){
           start_index = CRandom::GetRandomInteger( 0, dedisp.m_MWADataCube.m_Timesteps );
           printf("Random time index = %d\n",start_index);
        }

        // insert FRB :
        dedisp.InjectFRB( dyna_spec, start_index, gValue );
     }else{
        printf("WARNING : no FRB injected, option -n %d specfied - verify it was intentional\n",gNumberOfFRBs);
     }
    
     dyna_spec.WriteFits( gOutFitsFile.c_str() );
     printf("Dynamic spectrum with injected FRB saved to file %s\n",gOutFitsFile.c_str() );
  }else{
     unsigned long seed = DEFAULT_SEED + get_dttm(); // for testing keep the same + get_dttm();
     CRandom::Initialize( seed );
     CRandom rand;
     printf("Inserting %d FRBs using random seed = %d\n",gNumberOfFRBs,int(seed));
     printf("Using list of %d DMs and %d fluences\n",int(gDM_list.size()),int(gFluenceList.size()));

     FILE* out_f = fopen("generated.txt","w");
     fprintf(out_f,"# FITS DM Fluence[Jy ms] Flux[Jy] START-TIMEINDEX END-TIMEINDEX\n");
  
     for(int n=0;n<gNumberOfFRBs;n++){
        printf("--------------------------------------------- FRB %d ---------------------------------------------\n",n);
        dyna_spec.ReadFits( gInputFitsFile.c_str() );
        printf("\tRead read %s\n",gInputFitsFile.c_str() );
        
        int dm_index = rand.GetRandomInteger(0, gDM_list.size() );
        int fluence_index = rand.GetRandomInteger(0, gFluenceList.size() );
        int start_time_index = rand.GetRandomInteger(0, dyna_spec.GetXSize() );
        
        if( dm_index >= gDM_list.size() || fluence_index >= gFluenceList.size() ){
           printf("ERROR in code : index %d or %d larger than allowed !!!\n",dm_index,fluence_index);
           continue;
        }
        
        double dm = gDM_list[ dm_index ];
        double fluence = gFluenceList[ fluence_index ];
        double flux_jy = fluence / gIntTimeMiliseconds;
        char szPostfix[128];
        sprintf(szPostfix,"_genfrb%04d.fits",n);
        change_ext( gInputFitsFile.c_str(), szPostfix, gOutFitsFile, true );
        
        printf("\tDM                = %.2f [pc/cm^3] (index = %d)\n",dm,dm_index);
        printf("\tFluence           = %.2f [Jy ms]   (index = %d) -> flux density = %.4f [Jy]\n",fluence,fluence_index,flux_jy);
        printf("\tStart timde index = %d\n",start_time_index);
        printf("\tOutput fits file  = %s\n",gOutFitsFile.c_str());
        
        CDedispSearch::m_MinDM = dm;
        CDedispSearch::m_MaxDM = dm + 1; // due to change in CDedispSearch::CalcTemplatesExact( where double dm = m_MinDM + m_StepDM/2.00; while( dm <= m_MaxDM ){
                                
        dedisp.CalcTemplatesExact( dyna_spec );
        
        
        // insert FRB :
        // WARNING : this assumes no scattering (or the whole fluence contain in a 0.5sec (or other) time bin). 
        // If I am to add scattering here I will need to integrate flueance over all the 0.5sec bin (over which the signal is scattered in time)
        // and add these values of delta_flux_jy = fluence_bin_i / gIntTimeMiliseconds :
        int end_time_index=start_time_index;
        double injected_fluence = dedisp.InjectFRB( dyna_spec, start_time_index, flux_jy, false, gSNR, &end_time_index );
        if( gSNR > 0 ){
           // if injecting by SNR -> use the real value here:
           fluence = injected_fluence;
        }
                        
        
        dyna_spec.WriteFits( gOutFitsFile.c_str() );
        printf("\tDynamic spectrum with injected FRB saved to file %s\n",gOutFitsFile.c_str() );
        
        fprintf(out_f,"%s %.2f %.2f %.2f %d %d\n",gOutFitsFile.c_str(),dm,fluence,flux_jy,start_time_index,end_time_index); 
     }
     
     fclose( out_f );
  }

}
