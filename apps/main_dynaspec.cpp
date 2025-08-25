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

enum eTelescopeType {eUnknown=0, eMWA=1, eSKALOW_EDA2_AAVS2=2 };

eTelescopeType gTelescopeType = eMWA;
bool gChannelWidthSet=false;

string gOutputDir="dynamic_spectra/";
string gImageFileNameTemplate="wsclean_%d_timeindex%03d-%04d-I-dirty.fits";

int gCoarseChannels=24;
int gCoarseChannel=133;
double gStartFrequency=-1;
int gTimeSteps=-1;    // 592;
int gStartTimeIndex=0; 
double gTimeResolutionInSec=0.5;
int gObsID=1192530256;
double gThresholdInSigma=10;
int gBorder=20;
mystring gSubDirTemplate="wsclean_timeindex%03d";

bool gDumpSinglePixel = false;
int gBorderStartX=-1;
int gBorderStartY=-1;

int gBorderEndX=-1;
int gBorderEndY=-1;

// treating of zeros :
bool gZeros2NaNs = true;
double gZeroDistanceThreshold = 1e-9;

extern int gDebugLocal;

// order of channels / time directory structure :
enum eInputFitsFilesTypes  { eChannelDirFirst=0, eDefaultMWAFiles=1, eBlinkImager=2 };
// bool gChannelDirFirst=false;
eInputFitsFilesTypes gInputFitsFilesTypes = eDefaultMWAFiles;

void usage()
{
   printf("create_dynaspec -t TIME_STEPS -w (x_start,y_start)-(x_end,y_end) -p (X,Y)\n");
   printf("Description : program reads multiple images (in time and frequency) and saves resulting dynamic spectra to separate FITS files (for each image pixel)\n");
//    -n THRESHOLD_IN_SIGMA -b BORDER -a ENABLE_ALGO_TYPE -d OUTDIR -o OBSID -f FITS_FILENAME_TEMPLATE -r MAX_ALLOWED_RMS -D\n");   
   printf("\t-P option for reading PaCER BLINK images as they are in the format : start_time_1508442495_int_49_coarse_131_fine_ch00_image_real.fits\n");
   printf("\t-S START_TIME_INDEX : default = %d\n",gStartTimeIndex);   
   printf("\t-t TIME_STEPS : set <=0 to analyse all available timesteps\n");
   printf("\t-w (x_start,y_start)-(x_end,y_end) - do dump dynamic spectra of all pixels in this window\n");
   printf("\t-p (X,Y) - to dump a single pixel [default center pixel]\n");

//   printf("\t-a ENABLE_ALGO_TYPE : 1 - brute force, 2 - dynamic spectrum based also, 3 - both\n");
   printf("\t-o OBSID : obsid to be analysed [ default %d ]\n",gObsID);
   printf("\t-c CUT_THRESHOLD_IN_SIGMAS : threshold to cut values in dynamic spectrum in the dynaspec search [default %.2f]\n",CDedispSearch::m_ThresholdToCutInSigmas);   
   printf("\t-f FITS_FILENAME_TEMPLATE  : must contain %%d_timeindex%%03d-%%04d in the same order [default %s]\n",gImageFileNameTemplate.c_str());
   printf("\t-v : increase debug level by 1 [default debug level = %d]\n",gDebugLocal);
   printf("\t-r MAX_ALLOWED_RMS : default = %.4f\n",CDedispSearch::m_MaxRMSOnSingle);
   printf("\t-R EXCLUDE_REF_SOURCES : exclude reference sources [default %d]\n",CDedispSearch::m_bRejectRefSources);
   printf("\t-T Subdir template [default %s]\n",gSubDirTemplate.c_str());
   printf("\t-x Skip N first and last images [default %d]\n",CDedispSearch::m_SkipTimeSteps);
   printf("\t-D : enables de-dispersion debugging -> saves FITS files of pre-dedispersion images %.3fsec/1.28MHz [default %d]\n",gTimeResolutionInSec,CDedispSearch::m_bDebugDedispersion);
   printf("\t-C FIRST_COARSE_CHANNEL : first coarse channel [default %d]",gCoarseChannel);
   printf("\t-X TIME RESOLUTION IN SECONDS [default %.3f sec]\n",gTimeResolutionInSec);
   printf("\t-F : directory structure is CHANNEL/TIME/ (rather than default TIME/CHANNELS )\n");
   printf("\t-N N_CHAN : number of frequency channels [default %d]\n",gCoarseChannels);
   printf("\t-A BW_FINE_CH : frequency_resolution of a single channel in MHz [default %.2f MHz for MWA], for EDA2 = (32/27)*(400/512)/n_chan\n",MWA_COARSE_CHANNEL);
   printf("\t-B START_FREQ_MHZ : for MWA : coarse_channel*1.28-0.64, for EDA2/AAVS2 single channel : [ coarse_channel-16/27 ]*(400/512)\n");
   printf("\t-I TELESCOPE_TYPE : 1-MWA [default], 2 - EDA2/AAVS2, 0 - unknown [default %d]\n",gTelescopeType);
   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ho:l:m:s:t:n:b:a:d:c:f:vr:R:S:T:x:D:w:p:C:X:FN:A:B:I:P";
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

          case 'A':
            if( optarg ){
               MWA_COARSE_CHANNEL = atof( optarg );
               gChannelWidthSet = true;
            }
            break;

          case 'B':
            if( optarg ){
               gStartFrequency = atof( optarg );
            }
            break;

          case 'I':
            if( optarg ){
               gTelescopeType = (eTelescopeType)(atol( optarg ));
            }
            break;

          case 'd':
            if( optarg ){
               gOutputDir = optarg;
            }
            break;

          case 'C':
            if( optarg ){
               gCoarseChannel = atol( optarg );
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

          case 'F':
            // gChannelDirFirst = true;
            gInputFitsFilesTypes = eChannelDirFirst;
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

         case 'N':
            if( optarg ){
              gCoarseChannels = atol( optarg );
              MWA_COARSE_CHANNELS = gCoarseChannels;
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

         case 'X':
            if( optarg ){
               gTimeResolutionInSec = atof( optarg );
            }
            break;

         case 'w':
            if( optarg ){
               if( sscanf( optarg,"(%d,%d)-(%d,%d)",&gBorderStartX,&gBorderStartY,&gBorderEndX,&gBorderEndY )==4 ){
                  printf("Window correctly read , will save dynamic spectra of pixels in window (%d,%d) - (%d,%d)\n",gBorderStartX,gBorderStartY,gBorderEndX,gBorderEndY);
                  gDumpSinglePixel = false;
               }else{
                  printf("window %s is not correctly defined, exiting now\n",optarg);
                  exit(0);
               }
            }
            break;

         case 'p':
            if( optarg ){
               if( sscanf( optarg,"(%d,%d)",&gBorderStartX,&gBorderStartY )==2 ){
                  printf("Window correctly read , will dump dynamic spectrum of pixel (%d,%d)\n",gBorderStartX,gBorderStartY);
                  gDumpSinglePixel = true;
                  gBorderEndX = gBorderStartX;
                  gBorderEndY = gBorderStartY;
               }else{
                  printf("window %s is not correctly defined, exiting now\n",optarg);
                  exit(0);
               }
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
   
   printf("Creating output directory %s ...",gOutputDir.c_str() );
   MyFile::CreateDir( gOutputDir.c_str() );
   printf("OK\n");
      
   // checks :
   if( fabs(CDedispSearch::m_StepDM ) < 0.00000000001 ){
      printf("ERROR : too small DM step - the loop might become infinite -> specify larger value in -s option!!!!\n");
      exit(-1);            
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
}
 

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("Telescope type   = %d\n",gTelescopeType);
   printf("Start time index = %d\n",gStartTimeIndex);
   printf("Skip first images = %d\n",CDedispSearch::m_SkipTimeSteps);
   if( gDumpSinglePixel ){
      printf("Dump window = (%d,%d) - (%d,%d)\n",gBorderStartX,gBorderStartY,gBorderEndX,gBorderEndY);
   }else{
      printf("Dump single pixel = (%d,%d)\n",gBorderStartX,gBorderStartY);
   }
   printf("Input FITS files format = %d\n",int(gInputFitsFilesTypes));
   printf("\tOLD OPTION ;: Channel dir first = %d\n", (gInputFitsFilesTypes==eChannelDirFirst));
   printf("Number of freq. channels = %d\n",gCoarseChannels);
   printf("ObsID     = %d\n",gObsID);   

/*   printf("DM range  = %.2f - %.2f\n",CDedispSearch::m_MinDM,CDedispSearch::m_MaxDM);
   printf("DM step   = %.2f\n",CDedispSearch::m_StepDM);
   printf("Max allowed RMS = %.4f\n",CDedispSearch::m_MaxRMSOnSingle);
   printf("Border    = %d\n",gBorder);
   printf("Algorithms :\n");
   printf("\tBrute force = %d\n",CDedispSearch::m_bRunBruteForceAlgo);
   printf("\t\tParameters:\n");
   printf("\tDynaapsec   = %d\n",CDedispSearch::m_bRunDynSpecAlgo);
   printf("\t\tParameters:\n");
   printf("\t\tCutInSigma = %.2f sigmas\n",CDedispSearch::m_ThresholdToCutInSigmas);
   printf("Debug level = %d\n",gDebugLocal);
   printf("\t\tLog files used in de-dispersion = %d\n",CDedispSearch::m_bDebugDedispersion);
   printf("Exclude ref. sources = %d\n",CDedispSearch::m_bRejectRefSources);*/
   printf("Outdir    = %s\n",gOutputDir.c_str());
   printf("Timesteps = %d\n",gTimeSteps);
   printf("Subdir template = %s\n",gSubDirTemplate.c_str());
   if( gZeros2NaNs ){
      printf("Converting zeros ( < %e ) to NaNs\n",gZeroDistanceThreshold);
   }else{
      printf("WARNING : no conversion of zeros to NaNs\n");
   }
   printf("Start frequency      = %.6f [MHz]\n",gStartFrequency);
   printf("Start coarse channel = %d\n",gCoarseChannel);   
   printf("Time resolution      = %.6f [sec]\n",gTimeResolutionInSec);
   printf("Frequency resolution = %.6f [MHz]\n",MWA_COARSE_CHANNEL);
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

//  CDedispSearch frb_search( gObsID, gCoarseChannels, gTimeSteps, gStartTimeIndex );
//  int n_images = frb_search.m_MWADataCube.Read( gSubDirTemplate.c_str() , gImageFileNameTemplate.c_str() );
//  int events_count = frb_search.Run( gThresholdInSigma, gBorder, 6 );

  int check_image = 10;
  CMWADataCube cube_first( gObsID, gCoarseChannels, check_image ); // was 1 but usually does not exist so 10 should be safe !
  int read_images = 0;
  
  if( gInputFitsFilesTypes == eChannelDirFirst ){
     if( gTimeSteps < 0 ){
        std::vector<string> channel_fits_list;
        char szFitsList[128];
        sprintf(szFitsList,"%05d/fits_list",0);
      
        gTimeSteps = bg_read_list( szFitsList , channel_fits_list );
        gTimeSteps = channel_fits_list.size();
        printf("INFO : number of timestamps automatically read from list file %s to be %d\n",szFitsList,gTimeSteps);
     }
  
     read_images = cube_first.ReadChanTime( gSubDirTemplate.c_str() , gImageFileNameTemplate.c_str() );
  }else{
     if( gInputFitsFilesTypes == eBlinkImager ){
        int n_seconds = int(round(gTimeSteps*gTimeResolutionInSec));
        read_images = cube_first.ReadBlinkImages( gImageFileNameTemplate.c_str(), gTimeResolutionInSec, n_seconds, gCoarseChannel, 24 );
     }else{
        read_images = cube_first.Read( gSubDirTemplate.c_str() , gImageFileNameTemplate.c_str() );
     }
  }
  
  if( read_images > 0 ){
     printf("Read check image (%d) corretly -> getting parameters\n",check_image);
     std::vector< std::vector< CMWAFits* > > dynaspec_map;  
     int x_size = -1;
     int y_size = -1;
     CMWAFits* pCheckImage = cube_first.GetImage( 0, check_image-1 );
     if( pCheckImage ){     
        x_size = (pCheckImage)->GetXSize();
        y_size = (pCheckImage)->GetYSize();
     }else{
        printf("ERROR : could not read check image (%d) -> cannot continue, exiting\n",check_image);
        exit(-1);
     }
     printf("Auto-detected image size = %d x %d\n",x_size,y_size);
     std::vector< CMWAFits* > dynaspec_line;
     dynaspec_line.assign( x_size, NULL );
     dynaspec_map.assign( y_size, dynaspec_line );

     // initialise dynamic spectra for each (X,Y) pixel :     
     for( int y=0; y < y_size; y++){
        for( int x=0; x < x_size; x++){
           dynaspec_map[y][x] = NULL;
           
           // gBorderStartX,gBorderStartY,gBorderEndX,gBorderEndY
           if( y>=gBorderStartY && y<=gBorderEndY && x>=gBorderStartX && x<=gBorderEndX ){
              char szName[64];
              sprintf(szName,"pixel_%04d_%04d",x,y);
              dynaspec_map[y][x] = new CMWAFits( szName, gTimeSteps, gCoarseChannels );
              
              // init with NaNs :
              (dynaspec_map[y][x])->SetNaN();
           }           
        }
     }

     // filling dynamic spectra : 
     // read FITS files and fill (X,Y) pixels for a given (time,channel) in dynamic spectra 
     double prev_uxtime = -1.00;
     for(int start_timeindex=0;start_timeindex<gTimeSteps;start_timeindex++){
        // read 24 coarse channels for a given timestep :
        CMWADataCube cube( gObsID, gCoarseChannels, 1, start_timeindex );
        
        int read_images = 0;
        if( gInputFitsFilesTypes == eChannelDirFirst ){
           read_images = cube.ReadChanTime( gSubDirTemplate.c_str() , gImageFileNameTemplate.c_str(), false );
        }else{       
           if( gInputFitsFilesTypes == eBlinkImager ){
              int n_seconds = int(round(gTimeSteps*gTimeResolutionInSec));
              read_images = cube_first.ReadBlinkImages( gImageFileNameTemplate.c_str(), gTimeResolutionInSec, n_seconds, gCoarseChannel, 24 );
           }else{
              read_images = cube.Read( gSubDirTemplate.c_str() , gImageFileNameTemplate.c_str(), false );
           }
        }
        
        CMWAFits* pImage = cube.GetImage( 0, 0 );
        double unix_time = -1;
        if( pImage ){
           unix_time = pImage->GetUnixTime();
           
           if( prev_uxtime > 0 ){
              double diff = unix_time - prev_uxtime;
              
              if( fabs(diff) > 0.6 ){
                 printf("ERROR in code : image at start_timeindex = %d has unix_time = %.4f whilst previous unix_time = %.4f -> verify before using resulting dynamic spectra !!!\n",start_timeindex,unix_time,prev_uxtime);               
              }else{
                 printf("OK : start_timeindex = %d has unix_time = %.4f , previous unix_time = %.4f -> diff = %.4f sec\n",start_timeindex,unix_time,prev_uxtime,diff);
              }
           }        
           prev_uxtime = unix_time;
        }else{
           printf("WARNING : no image at start_timeindex = %d -> adding %.3f second to prev_uxtime = %.4f (if prev_uxtime>0) -> skipping\n",start_timeindex,gTimeResolutionInSec,prev_uxtime);
           if( prev_uxtime > 0 ){
              prev_uxtime = prev_uxtime + gTimeResolutionInSec;
              printf("WARNING : prev_uxtime := %.4f\n",prev_uxtime);
           }else{
              printf("WARNING : skipped start_timeindex = %d as prev_uxtime = %.2f (no images in the very beginning)\n",start_timeindex,prev_uxtime);
           }
           
           continue;
        }

        // set values of FITS image for all channels and given timestep :     
        for(int channel = 0; channel < gCoarseChannels ; channel++ ){
            for( int y=gBorderStartY; y <= gBorderEndY ; y++){
               for( int x=gBorderStartX; x <= gBorderEndX; x++){
                  CMWAFits* pDynaSpec = (dynaspec_map[y][x]);
           
                  if( pDynaSpec ){
                     double value = (cube[channel][0])->getXY(x,y); // cube has just one image for given channel and timestep - no need to [channel][start_timeindex]
                     if( gZeros2NaNs ){
                        if( value == 0.00 || fabs(value) < gZeroDistanceThreshold ){ // if 0 or consistent with 0.00 - floats and doubles might be tricky in this respect 
                           value = 0.00/0.00; // WARNING : FP_NAN does not really do anything but set 0.00
                        }
                     }
                     
                     pDynaSpec->setXY( start_timeindex, channel, value ); 
                  }else{
                     printf("ERROR : dynamic spectrum at pixel (%d,%d) has not been allocated - skipped\n",x,y);
                  }
               }
           }
        } // end of of loop over channels 
     } // end of timestep loop
     
     // saving dynamic spectra :
     printf("DEBUG : y_size = %d\n",y_size);
     for( int y=0; y < y_size; y++){
        for( int x=0; x < x_size; x++){
            CMWAFits* pFits = dynaspec_map[y][x];
            
            if( pFits ){
               // save mean spectrum to txt file :
               CBgArray mean_spectrum, rms_spectrum;
               pFits->MeanLines( mean_spectrum , rms_spectrum );
               char szOutMeanSpectrum[1024];
               sprintf(szOutMeanSpectrum,"%s/%04d_%04d.mean_spectrum",gOutputDir.c_str(),x,y);
               mean_spectrum.SaveToFile( szOutMeanSpectrum, NULL, &rms_spectrum );
               printf("Saved mean spectrum to file %s\n",szOutMeanSpectrum);               
                                          
               char szFileName[1024];
               sprintf(szFileName,"%s/%04d_%04d.fits",gOutputDir.c_str(),x,y);
               
               // add header keywords :
               pFits->SetKeyword( "CTYPE2", "Frequency" );
               pFits->SetKeyword( "CUNIT2", "MHz" );
               pFits->SetKeywordFloat( "CRPIX2", 0.5 );
               pFits->SetKeywordFloat( "CDELT2", (float)MWA_COARSE_CHANNEL );
               pFits->SetKeywordFloat( "CRVAL2", (float)gStartFrequency ); // to be un-hardcoded !
               pFits->SetKeyword( "CTYPE1","Time" );               
               pFits->SetKeyword( "CUNIT1","sec");
               pFits->SetKeyword( "CRPIX1", 1 );
               pFits->SetKeywordFloat( "CDELT1", (float)gTimeResolutionInSec );
               pFits->SetKeyword( "CRVAL1", 0.00 ); // time starts from 0.00
               
               if( pFits->WriteFits( szFileName ) ){
                   printf("ERROR : could not write dynamic spectrum %s\n",szFileName);
               }else{
                   printf("OK : dynamic spectrum %s written ok\n",szFileName);
               }
               
               // save also after subtracting mean spectrum :
               pFits->SubtractColumn( mean_spectrum );
               sprintf(szFileName,"%s/%04d_%04d_subtr.fits",gOutputDir.c_str(),x,y);
               if( pFits->WriteFits( szFileName ) ){
                   printf("ERROR : could not write dynamic spectrum %s\n",szFileName);
               }else{
                   printf("OK : dynamic spectrum %s written ok\n",szFileName);
               }
            }
        }
     }
  }else{
     printf("ERROR : could not read the test image\n");
     exit(-1);
  }

  
//  CMWAFits dynspec( "dynaspec", m_MWADataCube.m_Timesteps, m_MWADataCube.m_Channels );
    
//  CMWADataCube data_cube( gObsID, gCoarseChannels, gTimeSteps );
//  int n_images = data_cube.Read();


   time_t end_time = get_dttm();
   printf("FRB finder finished at unixtime = %d, took %d seconds\n",(int)end_time,(int)(end_time-start_time));
}
