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

enum eTelescopeType {eUnknown=0, eMWA=1, eSKALOW_EDA2_AAVS2=2 };
eTelescopeType gTelescopeType = eMWA;
bool gChannelWidthSet=false;
double gStartFrequency=-1;


int gCoarseChannels=24;
int gCoarseChannel=133;
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
   printf("dynspec_search FITS_FILE_DYNA-SPEC OUT_FITS_FILE -l MIN_DM -m MAX_DM -s DM_STEP -t TIME_STEPS -n THRESHOLD_IN_SIGMA -b BORDER -o OBSID -a AVERAGE_N_TIMESTEPS\n");
   printf("-a AVERAGE_N_TIMESTEPS : to enable averaging of N timesteps before de-dispersion [default %d (<=1 is disabled)]\n",gAverageNTimeSteps);
//   printf("-c Threshold_on_SNR : threshold for SNR [default = %.2f]\n",CDedispSearch::m_SNRThreshold);
//   printf("-x Threshold_for_cleaning : threshold to ignore pixels with lower value [default %.2f]\n",CDedispSearch::m_SNRThreshold);
   printf("-n THRESHOLD_IN_SIGMA : threshold to cut pixels in the dynamic spectrum [default %.2f]\n",gThresholdInSigma);
   printf("-N THRESHOLD_IN_SIGMA for SNR image :  threshold to cut pixels in the dynamic spectrum [default %.2f]\n",gThresholdInSigmaForSNR);
   printf("-o OBSID : observation ID [default %d]\n",gOBSID);
   printf("-d DEBUG_LEVEL\n");
   printf("-L LOG_BELOW_DM [default = %e = turned of]\n",gLogBelowDM);
   printf("-X : exclude maximum value in dispersion sweep to ignore RFI polluted channels [default disabled]\n");
   printf("-w 0/1 : turn off or on weigtiing of pixels by (1.0/RMS) when de-dispersing [default %d]\n",gDoWeigthByRMS);
   printf("-r : use old rms map (more conservative and usually higher than the new one [default]\n");
   printf("-A : auto metafits fits\n");
   printf("-D DEBUG_DM\n");
   printf("-T : use global threshold (captial T for global, default is local RMS per channel)\n");
   printf("-x : transposed FITS files for dynamic spectrum is horizontal time on X-axis and frequency on Y-axis\n");
   printf("-F START-END : range of channels to be flagged, for example -F 0-10\n");
   printf("-f : enable auto-flagging\n");
   printf("-S SEARCH_TYPE : DTS search type 0 - no search , 1 - standard (OLD) DTS search [default]\n");
   printf("-P 0/1 : allow for dispersion sweeps starting before the start of dynamic spectrum (1) , 0 - to not allow for this case [default %d]\n",CDedispSearch::m_bAllowPreStartArrivals);
   printf("\t-C BW_FINE_CH : frequency_resolution of a single channel in MHz [default %.2f MHz for MWA], for EDA2 = (32/27)*(400/512)/n_chan\n",MWA_COARSE_CHANNEL);
   printf("\t-B START_FREQ_MHZ : center of the first channel (can be taken from FITS header) for MWA : coarse_channel*1.284, for EDA2/AAVS2 single channel : [ coarse_channel-16/27 ]*(400/512)\n");
   printf("\t-E COARSE CHANNEL [default %d]\n",gCoarseChannel);
   printf("\t-I TELESCOPE_TYPE : 1-MWA [default], 2 - EDA2/AAVS2, 0 - unknown [default %d]\n",gTelescopeType);

   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "a:ho:l:m:s:t:n:N:b:d:o:L:Xw:rAD:TxF:fS:C:B:I:E:P:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'a':
            if( optarg ){
               gAverageNTimeSteps = atol( optarg );
            }
            break;

         case 'A':
            gAutoMetafits = true;
            break;

         case 'C':
           if( optarg ){
               MWA_COARSE_CHANNEL = atof( optarg );
               gChannelWidthSet = true;
            }
            break;

         case 'I':
           if( optarg ){
              gTelescopeType = (eTelescopeType)(atol( optarg ));
           }
           break;

         case 'E':
           if( optarg ){
              gCoarseChannel = (atol( optarg ));
           }
           break;

         case 'B':
            if( optarg ){
               gStartFrequency = atof( optarg );
            }
            break;

         case 'D':
            CDedispSearch::m_DebugDM = atof( optarg );
            break;

         case 'o':
            if( optarg ){
               gOBSID = atol( optarg );
            }
            break;

/*         case 'c':
            if( optarg ){
               CDedispSearch::m_SNRThreshold  = atof( optarg );
            }
            break;*/

         case 'n':
            if( optarg ){
               gThresholdInSigma = atof( optarg );
               CDedispSearch::m_SNRThreshold = gThresholdInSigma;
            }
            break;

         case 'N':
            if( optarg ){
               gThresholdInSigmaForSNR = atof( optarg );
               // CDedispSearch::m_SNRThreshold = gThresholdInSigma;
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

         case 'r':
            gDivideByMapNew = false;
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

         case 'f':
            gAutoFlagging = true;
            break;

         case 'F':
            if( optarg ){
                gFlaggedRanges.push_back( string(optarg) );
            }
            break;

         case 'L':
            if( optarg ){
              gLogBelowDM = atof( optarg );
            }
            break;

         case 'P':
            if( optarg ){
              CDedispSearch::m_bAllowPreStartArrivals = ( atol( optarg ) > 0 );
            }
            break;

         case 'S':
            if( optarg ){
              gSearchType = (eDTSSearchType)( atol( optarg ) );
            }
            break;

         case 'T':
            CMWAFits::m_bUseLocalThreshold = false;
            break;

         case 'X':
            gIgnoreMaxValue = true;
            break;

         case 'w':
            gDoWeigthByRMS = ( atol( optarg ) > 0 );
            break;

         case 'x':
            gTransposedFits = true;
            break;

         case 'h':
            usage();
            break;
         default:  
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
   
   gBadChannelRanges.clear();
   CFlagRanges tmp;
   for(int i=0;i<gFlaggedRanges.size();i++){
      if( sscanf(gFlaggedRanges[i].c_str(),"%d-%d",&(tmp.m_Start),&(tmp.m_End)) == 2 ){
         gBadChannelRanges.push_back( tmp );
         printf("PARAMETERS : flagged channel range %d - %d added to list of flagged channels\n",tmp.m_Start,tmp.m_End);
      }else{
         printf("ERROR : could not parse flagged channel range %s -> ignored\n",gFlaggedRanges[i].c_str());         
      }
   }

   if( gStartFrequency < 0 ){
      if( gTelescopeType == eSKALOW_EDA2_AAVS2 ){
         // lower edge of the band - start of the very first fine channel :
         // gStartFrequency = ( double(gCoarseChannel) - 16.00/27.00 )*(400.00/512.00);
         // center of the very first fine channel - it was the same previously :
         if( !gChannelWidthSet ){
            MWA_COARSE_CHANNEL = (32.00/27.0)*(400.00/512.00) / gCoarseChannels;
            printf("EDA2/AAVS2 : fine channel width automatically calculated to be %.6f MHz\n",MWA_COARSE_CHANNEL);
         }
         gStartFrequency = double(gCoarseChannel)*(400.0/512.0) - (32.00/27.0)*(400.00/512.00)*(1.00/2.00) + MWA_COARSE_CHANNEL/2.00;
         printf("EDA2/AAVS2 : start frequency automatically calculated to be %.6f MHz\n",gStartFrequency);
      }else{
         gStartFrequency = (gCoarseChannel*1.28);
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
   printf("Telescope type   = %d\n",gTelescopeType);
   printf("Input fits file  = %s (transposed = %d)\n",gInputFitsFile.c_str(),gTransposedFits);
   printf("Output fits file = %s\n",gOutFitsFile.c_str()); 
   printf("Output reg file  = %s\n",gOutCandidatesFile.c_str());
   printf("OBSID     = %d\n",gOBSID);  
   printf("DM range  = %.2f - %.2f\n",CDedispSearch::m_MinDM,CDedispSearch::m_MaxDM);
   printf("DM step   = %.2f\n",CDedispSearch::m_StepDM);
   printf("Timesteps = %d\n",gTimeSteps);
   printf("Border    = %d\n",gBorder);
   printf("Threshold on SNR   = %.2f (or %.2f)\n",gThresholdInSigmaForSNR,CDedispSearch::m_SNRThreshold);
   printf("Threshold to clean = %.2f\n",gThresholdInSigma);
   printf("Low below DM = %e\n",gLogBelowDM);
   printf("Ignore max value = %d\n",gIgnoreMaxValue);
   printf("Weithing enabled = %d\n",gDoWeigthByRMS);
   printf("Average N timesteps = %d\n",gAverageNTimeSteps);
   printf("Use local threshold in DTS search = %d\n",CMWAFits::m_bUseLocalThreshold);   
   printf("Do search = %d\n",gSearchType);   

   printf("Coarse channel       = %d\n",gCoarseChannel);   
   printf("Start frequency      = %.20f [MHz]\n",gStartFrequency);
//   printf("Start coarse channel = %d\n",gCoarseChannel);   
//   printf("Time resolution      = %.6f [sec]\n",gTimeResolutionInSec);
   printf("Frequency resolution = %.20f [MHz]\n",MWA_COARSE_CHANNEL);
   printf("Allow for pre-start arrival times = %d\n",CDedispSearch::m_bAllowPreStartArrivals);  
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

  parse_cmdline(argc-2,argv+2);
  print_parameters();


  // const char* change_ext(const char* name,const char* new_ext,string& out);  
  string gOutRmsMapFile, gOutSnrSeriesFile, gOutSnrSeriesRegFile, gOutRmsMapNewFile;
  change_ext( gOutFitsFile.c_str(), "reg", gOutCandidatesFile );
  change_ext( gOutFitsFile.c_str(), "_count.fits", gOutCountFitsFile, true );
  change_ext( gOutFitsFile.c_str(), "_rmsmap.fits", gOutRmsMapFile, true );
  change_ext( gOutFitsFile.c_str(), "_rmsmapnew.fits", gOutRmsMapNewFile, true );
  change_ext( gOutFitsFile.c_str(), "_snr.fits", gOutSnrSeriesFile, true );
  change_ext( gOutFitsFile.c_str(), "_snr.reg", gOutSnrSeriesRegFile, true );
  
  char szTmp[64];
  sprintf(szTmp,"_avg%02dtimesteps.fits",gAverageNTimeSteps);
  change_ext( gOutFitsFile.c_str(), szTmp, gOutAvgDynSpec, true );


  CMWAFits dyna_spec( gInputFitsFile.c_str() );
  dyna_spec.ReadFits( gInputFitsFile.c_str(), 0, 1, 0, gTransposedFits );

  // TODO : add parameters and turn on/off 
  if( gAutoFlagging ){
     printf("INFO : automatically flagging RFI affected channels\n");
     dyna_spec.AutoFlagChannels(0.5); // 0.5sigma threshold
  }else{
     printf("WARNING : no automatic RFI flagging requested\n");
     // dyna_spec.m_IsFlagged.assign( dyna_spec.GetYSize(), false );
     dyna_spec.InitFlagsByRanges( gBadChannelRanges );
  }
  
  if( gAverageNTimeSteps > 1 ){
     dyna_spec.avg_in_time(  gAverageNTimeSteps , gOutAvgDynSpec.c_str() );
  }

  double delta_freq=dyna_spec.delta_freq;
  if( gChannelWidthSet ){
     delta_freq = MWA_COARSE_CHANNEL;
  }
  
  double start_freq=dyna_spec.start_freq;
  if( gStartFrequency > 0 ){
     start_freq = gStartFrequency;
  }


  CDedispSearch dedisp( gOBSID, dyna_spec.GetYSize(), dyna_spec.GetXSize() );
//  dedisp.m_MWADataCube.m_StartUnixTime = 
  dedisp.m_MWADataCube.m_Timesteps = dyna_spec.GetXSize();
  dedisp.m_MWADataCube.m_Channels  = dyna_spec.GetYSize();
  dedisp.m_MWADataCube.m_FreqUpperMHz = start_freq + dyna_spec.GetYSize()*delta_freq - delta_freq/2.00; // upper end of the band (edge of highest channel - NOT CENTRE !)
  printf("DEBUG : STARTFREQ = %.6f MHz, ENDFREQ = %.6f MHz\n",start_freq,dedisp.m_MWADataCube.m_FreqUpperMHz);

  if( gAutoMetafits ){
     printf("DEBUG : reading of metafits file is not required\n");
     dedisp.m_MWADataCube.GenerateMetaData( dyna_spec, dedisp, gOBSID, start_freq, delta_freq );
  }else{
     dedisp.m_MWADataCube.ReadMetaData( gOBSID );
  }

  std::vector<double> freq_list;
  if( dedisp.m_MWADataCube.get_freq_list( freq_list ) != dedisp.m_MWADataCube.m_CoarseChannels.size() ||  dedisp.m_MWADataCube.m_CoarseChannels.size() <= 0 ){
     printf("ERROR : could not create the list of observing frequencies\n");
     exit(-1);
  }


  double inttime = dyna_spec.inttime;
  if( dedisp.m_MWADataCube.m_Metafits ){
     inttime = (dedisp.m_MWADataCube.m_Metafits)->inttime;
  }
  dedisp.CalcTemplatesExact( dyna_spec );
  
  
  int min_good_pixels = 12; // 20;
//  dedisp.AnalyseDynSpec( dyna_spec, 0, 0, gThresholdInSigma, &candidates_file, 1, min_good_pixels, &freq_list );

  CMWAFits dedispersed_series(gOutFitsFile.c_str()), count_map, rms_map_new;
/*
  int GetDedispersedSeries( CMWAFits& dynspec, CMWAFits& dedispersed_series, CMWAFits& out_count_map, std::vector<double>* p_channels_list=NULL, double threshold_in_sigma=5.00, 
                             MyOFile* candidates_file=NULL, int max_radius=1, double LogBelowDM=-1e6 );   
*/    
  // rms_map_new is map of RMS-s calculated as Sum of RMS_channel^2 for all pixels in the DM-sweep
  int max_count = dedisp.GetDedispersedSeries( dyna_spec, dedispersed_series, count_map, rms_map_new, NULL, 5.00, NULL, 1, gLogBelowDM, gIgnoreMaxValue, gDoWeigthByRMS ); //  std::vector<double>* p_channels_list, double threshold_in_sigma, MyOFile* candidates_file, int max_radius )
  printf("INFO : m_MaxDispersiveSweep = %.6f [sec], inttime = %.6f sec , maximum number of pixels = %d\n",dedisp.m_MaxDispersiveSweep,inttime,max_count);
  dedispersed_series.WriteFits( gOutFitsFile.c_str() );
  count_map.WriteFits( gOutCountFitsFile.c_str() );
  
  // adding saving of MEAN and RMS per channel from DTS :
  CBgArray mean_channel_dts,rms_channel_dts;
  dedispersed_series.MeanLines( mean_channel_dts , rms_channel_dts );
  string szOutTmpFile;
  change_ext( dyna_spec.GetFileName(), "_rms_per_dm_DTS.txt", szOutTmpFile, true );
  rms_channel_dts.SaveToFile( szOutTmpFile.c_str() );
  change_ext( dyna_spec.GetFileName(), "_mean_per_dm_DTS.txt", szOutTmpFile, true );
  mean_channel_dts.SaveToFile( szOutTmpFile.c_str() );


  if( gSearchType == eNormalSearch ){
     // calculate what can be maximum dispersion delay in pixels from dedisp.m_MaxDispersiveSweep and time resolution:
     // int n_max_steps = dedisp.m_MaxDispersiveSweep / inttime; // dedispersed_series.cdelt1;
     // int max_count = sqrt( (n_max_steps)*(n_max_steps) + dyna_spec.GetYSize()*dyna_spec.GetYSize() );
     printf("DEBUG : max_disp_sweep = %.6f [sec] -> max_count = %d\n",dedisp.m_MaxDispersiveSweep,max_count);
     
  
     // RMS per count is calculated as RMS value corresponding to number of pixels in DM-sweep
     CBgArray rms_per_count;  
     dedispersed_series.FindSources( gThresholdInSigma, 0, gOutCandidatesFile.c_str(), true, 60, -1, &count_map, &rms_per_count, max_count ); // WARNING - HARDCODING !!! // 2021-04-14 : was 17200 -> -1 (to analyse all timestamps)
  
     // calculate RMS map using RMS per count lookup table :
     CMWAFits rms_map( "rms_map.fits", dedispersed_series.GetXSize(), dedispersed_series.GetYSize() );  
     for(int y=0;y<rms_map.GetYSize();y++){
        for(int x=0;x<rms_map.GetXSize();x++){
           int count = count_map.getXY( x, y );
           if( count < rms_per_count.size() ){
              double rms_value = rms_per_count[count];
              rms_map.setXY( x, y, rms_value );
           }else{
              printf("ERROR : rms value not calculated for count = %d -> setting NaN at (%d,%d) in rms_map.fits\n",count,x,y);
              rms_map.setXY( x, y, (0.00/0.00) );
           }
        }
     }  
  
     // calculated from group of paths with the same number of pixels 
     rms_map.WriteFits( gOutRmsMapFile.c_str() );  
     printf("RMS map saved to file %s\n",gOutRmsMapFile.c_str() );
  
     // calculated as a weighted sum of RMS^2 values where RMS is calculated as a function of frequency ...
     // TODO : it is very underestimated - implement Jun's idea :
     //    1/ select 100 pixels or so
     //    2/ calculate DTS for each of them
     //    3/ For each pixel in DTS calculate RMS out of these 100 DTS images at given (x,y) position 
     //    4/ TEST on Gemma's images from the paper as a double/cross-check !  
     rms_map_new.WriteFits( gOutRmsMapNewFile.c_str() );
     printf("RMS map-NEW saved to file %s\n",gOutRmsMapNewFile.c_str() );

     //  if( gDivideByMapNew ){
     if( false ){ // using RMS-MAP new has been disabled because it is very much underestimating RMS -> too many high SNR event candidates !
        printf("INFO : dividing by new rms_map - should be ok now\n");
        dedispersed_series.Divide( rms_map_new );
     }else{  
        printf("WARNING : dividing by old type of rms_map - please verify if OK\n");     
        dedispersed_series.Divide( rms_map );    
     }
     dedispersed_series.WriteFits( gOutSnrSeriesFile.c_str() );
     printf("Saved SNR FITS file %s\n",gOutSnrSeriesFile.c_str() );  

     // int FindSourcesSNR( double threshold_snr=5.00, int border = 20, const char* szOutRegFile="sources_snr.reg", int minX=0, int maxX=-1 )  
     // gOutSnrSeriesRegFile
     // find time range :
     double max_time = dyna_spec.GetXSize()*dyna_spec.inttime;
     printf("INFO : max_time = %.4f [sec] , max_dispersive_sweep = %.4f [sec]\n",max_time,dedisp.m_MaxDispersiveSweep);
     int start_time_index=-1,end_time_index=dedispersed_series.GetXSize();
     for(int t=0;t<dedispersed_series.GetXSize();t++){
        double time = dedispersed_series.m_StartTime + t*dedispersed_series.inttime;
        if( time >= 0 && start_time_index<0 ){
           start_time_index = t;
           printf("\tSet start time index = %d\n",start_time_index);
        }
        if( time >= (max_time-dedisp.m_MaxDispersiveSweep) ){
           end_time_index = t;
           printf("\tSet end time index = %d\n",end_time_index);
           break;
        }
     }
     printf("INFO : checking range of times %d - %d\n",start_time_index,end_time_index);
     dedispersed_series.FindSourcesSNR( gThresholdInSigmaForSNR, 5, gOutSnrSeriesRegFile.c_str(),  start_time_index+2, end_time_index-2, &count_map ); // 2021-04-14 : was 17200 -> -1 (to analyse all timestamps)
   }else{
      printf("WARNING : searach for FRBs in DTS is not required\n");      
   }
}
