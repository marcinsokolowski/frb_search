#include "dedisp_search.h"
#include <math.h>
#include <myfile.h>
#include <mystring.h>
#include "myprogress.h"

// test only 
// #include <myfits.h>

#define DEDISPERSION_VERSION "20191102"

// #include <bits/stdc++.h> 


bool   CDedispSearch::m_bAllowPreStartArrivals = true; // calculate DTS with negative arrival times (max negative value = max(dipersive delay))
bool   CDedispSearch::m_bAllowNegativeDM = true;
double CDedispSearch::m_MinDM  = 200;
double CDedispSearch::m_MaxDM  = 1000;
double CDedispSearch::m_StepDM = 1;
double CDedispSearch::m_SNRThreshold = 5;
double CDedispSearch::m_ThresholdToCutInSigmas = -1000; // threshold to cut values below this threshold from the dynamic spectrum in 
                                                        // the FindDynaSpec algorithm 

double CDedispSearch::m_MinRMSOnSingle  = 0.00001;
double CDedispSearch::m_MaxRMSOnSingle  = 5.00; // maximum allowed RMS on a single 0.5 sec (or other) image 
int    CDedispSearch::m_MinUsedChannels = 8;   // minimum number of used channels , was 10, but sqrt(10) = 3.16 and sqrt(8) = 2.82, so it's not such a huge difference in noise suppression so 
                                               // I can try to allow only 8 channels 
int    CDedispSearch::m_MinUsedImagesInDedisp = 24; // use 24 coarse channel number as a minimum to search for transients
int    CDedispSearch::m_SkipTimeSteps   = 0;    // how many timesteps to skip , was 10 and 8, but now I will relay on automatic detection of ZERO images !

string CDedispSearch::m_ResultsPath = "candidates/";
int CDedispSearch::m_bRunBruteForceAlgo=1;
int CDedispSearch::m_bRunDynSpecAlgo=0;

int CDedispSearch::m_bRejectRefSources = 1;
bool CDedispSearch::m_bForceOverwrite = false;

// DEBUG :
bool CDedispSearch::m_bDebugDedispersion = false;
double CDedispSearch::m_DebugDM = -1;

int CDedispSearch::find_value_double( vector<double>& tab, double value, double error /*=0.000001*/ )
{
        int pos=0;
        for( vector<double>::iterator i=tab.begin();i!=tab.end();i++,pos++){
                if( fabs( (*i)-value ) <= error ){
                   return pos;
                }
        }
        return -1;
}



CDispersedPixels::CDispersedPixels()
{}

void CDispersedPixels::save_text_file( const char* filename )
{
   FILE* outf = fopen( filename , "w" );
   for(int i=0;i<size();i++){
      CPixel& pix = (*this)[i];      
      fprintf(outf,"%.8f %.8f %.8f\n",pix.exact_uxtime,pix.exact_freq_mhz,pix.flux);
   }
   fclose( outf );
}

void CDispersedPixels::save_reg_file( const char* filename )
{
   FILE* outf = fopen( filename , "w" );
   for(int i=0;i<size();i++){
      CPixel& pix = (*this)[i];      
      fprintf(outf,"circle %d %d 1 # color=green\n",pix.time,pix.cc);
   }
   fclose( outf );
}

bool CDispersedPixels::add_unique( CPixel& new_pixel )
{
   bool bFound=false;
   for(int i=0;i<size();i++){
      if( (*this)[i].time == new_pixel.time && (*this)[i].cc == new_pixel.cc ){
         // found not addning
         return false;
      }
   }
   
   push_back( new_pixel );
   return true;
}


void CEventList::add( int start_timeidx, int x, int y, eAlgorithmType_T _algo_type, double _dm, double _snr, double _flux )
{
   CEvent new_event( start_timeidx, x, y, _algo_type, _dm, _snr, _flux );
   push_back( new_event );
}

double CDedispSearch::dispersion_milisec( double freq1_mhz, double freq2_mhz, double dm )
{
   double freq1_ghz = freq1_mhz / 1000.0;
   double freq2_ghz = freq2_mhz / 1000.0;

   double delta_t_ms = 4.15 * dm * ( 1.00/(freq2_ghz*freq2_ghz) - 1.00/(freq1_ghz*freq1_ghz) );
   return delta_t_ms;
}   

double CDedispSearch::dispersion_freq_mhz( double freq1_mhz, double time_value_ms, double dm )
{
   double freq_mhz = freq1_mhz;
   
   if( dm > 0.0000001 ){ // really >0 :
      double freq1_ghz = freq1_mhz / 1000.00;
      double freq_ghz = 1.00 / sqrt( 1.00/(freq1_ghz*freq1_ghz) + time_value_ms / (dm*4.15) ); // correct formula other way around    
      freq_mhz = freq_ghz * 1000.00;  
   }
   
   return freq_mhz;
}   

CDedispSearch::CDedispSearch( int obsid, int n_channels, int timesteps, int start_timeindex )
: m_MWADataCube( obsid, n_channels, timesteps, start_timeindex ), m_DedispersedImagesIndex(0), m_MaxDispersiveSweep(0.00)
{
}

CDedispSearch::~CDedispSearch()
{
   CleanDedispersed(1);
}   

void CDedispSearch::CleanDedispersed( int bDelete )
{
   for(int i=0;i<m_DedispersedImages.size();i++){
       if( bDelete > 0 ){
           delete m_DedispersedImages[i];
       }else{
          (m_DedispersedImages[i])->SetValue(0);
       }
   }
   m_DedispersedImagesIndex = 0;
   // m_DedispersedImages.clear();  
}

void CDedispSearch::CalcReferenceImages()
{
    // de-disperse at 0 DM to caculate mean image :
    // Dedisperse( 0.00 );
 
    // clean median and mean images if required :
    if( m_MWADataCube.m_pMeanImage ){
       (m_MWADataCube.m_pMeanImage)->SetValue( 0.00 );
       (m_MWADataCube.m_pMeanImage)->reset_lines_counter();
    }else{
       printf("ERROR : m_MWADataCube.m_pMeanImage = NULL -> cannot continue\n");
       exit(-1);
    }

    if( m_MWADataCube.m_pRMSImage ){
       (m_MWADataCube.m_pRMSImage)->SetValue( 0.00 );
       (m_MWADataCube.m_pRMSImage)->reset_lines_counter();
    }else{
       printf("ERROR : m_MWADataCube.m_pRMSImage = NULL -> cannot continue\n");
       exit(-1);
    }
    
    for(int timeindex=0;timeindex<m_MWADataCube.m_Timesteps;timeindex++){
        // double unixtime = m_MWADataCube.m_StartUnixTime + start_timeindex * inttime;

        for(int cc=0;cc<m_MWADataCube.m_Channels;cc++){
            CMWAFits* pImage = m_MWADataCube.GetImage( cc, timeindex );
            
            if( pImage ){
                pImage->GetStat( CMWAFits::m_Border );
                
                printf("CDedispSearch::CalcReferenceImages : %s timeinex=%d, cc=%d  , RMS = %.4f\n",pImage->GetFileName(),timeindex,cc,pImage->m_RMS);
                
                if( pImage && pImage->m_ReadStatus == 0 && ( fabs(pImage->m_Mean)>0.0001 /* || fabs(pImage->m_RMS)>0.0001 */ || fabs(pImage->m_Max)>0.0001 ) ){
                    // add to mean image :
                    if ( pImage->m_RMS > m_MinRMSOnSingle && pImage->m_RMS < m_MaxRMSOnSingle ){
                       // only good images are included into MEAN/RMS :
                       (*m_MWADataCube.m_pMeanImage) += (*pImage);
                       (m_MWADataCube.m_pMeanImage)->inc_lines_counter();
                    
                       (m_MWADataCube.m_pRMSImage)->add_squared( *pImage );                    
                       (m_MWADataCube.m_pRMSImage)->inc_lines_counter();                    
                    
                        m_MWADataCube.m_ImageFlags[cc][timeindex] = 1;
                    }else{
                       if( gDebugLocal >= 0 ){printf("WARNING : timeinex=%d, cc=%d flagged and not used for MEAN/RMS calculation\n",timeindex,cc);}
                    }
                }else{
                   if( gDebugLocal >= 0 ){printf("WARNING : timeinex=%d, cc=%d skipped due to missing or zero image\n",timeindex,cc);}
                }
            }
        }
    }


    if( m_MWADataCube.m_pMedianImage ){
       (m_MWADataCube.m_pMedianImage)->SetValue( 0.00 );
    }else{
       printf("ERROR : m_MWADataCube.m_pMedianImage = NULL -> cannot continue\n");
       exit(-1);
    }

    if( m_MWADataCube.m_pMeanImage ){
        m_MWADataCube.m_pMeanImage->Divide( (m_MWADataCube.m_pMeanImage)->get_lines_counter() );
       
        m_MWADataCube.m_pMeanImage->SetKeyword( "N_IMAGES", (m_MWADataCube.m_pMeanImage)->get_lines_counter()  );
        
        char szOutFullPath[1024];
        sprintf(szOutFullPath,"%s/%s",m_ResultsPath.c_str(),m_MWADataCube.m_pMeanImage->GetFileName());
        m_MWADataCube.m_pMeanImage->SetFileName( szOutFullPath );
        if( m_MWADataCube.m_pMeanImage->WriteFits( szOutFullPath ) ){
            printf("ERROR : could not write mean image to file %s\n",szOutFullPath);
        }

        if( m_MWADataCube.m_pRMSImage ){
           sprintf(szOutFullPath,"%s/%s",m_ResultsPath.c_str(),m_MWADataCube.m_pRMSImage->GetFileName() );
           // m_MWADataCube.m_pRMSImage->Divide( (m_MWADataCube.m_pRMSImage)->get_lines_counter() );
           m_MWADataCube.m_pRMSImage->RMS( (m_MWADataCube.m_pRMSImage)->get_lines_counter(), (*m_MWADataCube.m_pMeanImage) );           
       
           m_MWADataCube.m_pRMSImage->SetKeyword( "N_IMAGES", (m_MWADataCube.m_pRMSImage)->get_lines_counter()  );
           m_MWADataCube.m_pRMSImage->SetFileName( szOutFullPath );
           if( m_MWADataCube.m_pRMSImage->WriteFits( szOutFullPath ) ){
               printf("ERROR : could not write mean image to file %s\n",m_MWADataCube.m_pRMSImage->GetFileName());
           }
        }
    }
       
    // TODO calculate median :
    if( m_MWADataCube.m_pMedianImage ){       
        printf("Calculating median of good / non-flagged images ...");fflush(stdout);
    
        double* sort_table = new double[m_MWADataCube.m_Timesteps*m_MWADataCube.m_Channels];
        int min_image_count = 1e6;
          
        for(int y=0;y<m_MWADataCube.m_pMedianImage->GetYSize();y++){
           printf("Progress : %d / %d\n",y,m_MWADataCube.m_pMedianImage->GetYSize());
           for(int x=0;x<m_MWADataCube.m_pMedianImage->GetXSize();x++){

              int n_image_count=0;
              for(int timeindex=0;timeindex<m_MWADataCube.m_Timesteps;timeindex++){
                  for(int cc=0;cc<m_MWADataCube.m_Channels;cc++){
                      if( m_MWADataCube.m_ImageFlags[cc][timeindex] > 0 ){
                          CBgFits* pImage = m_MWADataCube.GetImage( cc, timeindex );
                              
                          if( pImage ){
                              sort_table[n_image_count] = pImage->valXY(x,y);                          
                              n_image_count++;                                  
                          }
                      }
                  }
              }
              if( n_image_count > 0 ){
                  my_sort_float( sort_table, n_image_count );
                  double median = sort_table[ n_image_count / 2 ];
                  m_MWADataCube.m_pMedianImage->setXY( x, y, median );
                      
                  if( n_image_count < min_image_count ){
                     min_image_count = n_image_count;
                  }
              }
          }
      }
      delete [] sort_table;
       
      m_MWADataCube.m_pMedianImage->SetKeyword( "N_IMAGES", min_image_count );
      char szOutFullPath[1024];
      sprintf(szOutFullPath,"%s/%s",m_ResultsPath.c_str(),m_MWADataCube.m_pMedianImage->GetFileName() );
      m_MWADataCube.m_pMedianImage->SetFileName( szOutFullPath );
      if( m_MWADataCube.m_pMedianImage->WriteFits( szOutFullPath ) ){
          printf("ERROR : could not write mean image to file %s\n",m_MWADataCube.m_pMedianImage->GetFileName());
      }
      
      printf("Median done\n");fflush(stdout);
   }
   
   m_MWADataCube.FindReferenceSources();
}

// 20190909 - > changed to _OLD and _NEW to normal Dedisperse( double dm )
void CDedispSearch::Dedisperse_OLD( double dm )
{
    // prepare memory for de-dispersed images :
    if( m_DedispersedImages.size() < m_MWADataCube.m_Timesteps ){
        printf("DEBUG : getting image from timestep = %d , channel = %d\n",m_MWADataCube.m_FirstCorrectTimestep, m_MWADataCube.m_FirstCorrectChannel );
        CMWAFits* pImage = m_MWADataCube.GetImage( m_MWADataCube.m_FirstCorrectChannel , m_MWADataCube.m_FirstCorrectTimestep );
        
        if( pImage ){
            for(int start_timeindex=0;start_timeindex<m_MWADataCube.m_Timesteps;start_timeindex++){
                if( start_timeindex >= m_DedispersedImages.size() ){
                    char szName[512];
                    sprintf(szName,"starttimdeindex_%04d_dm%.2f.fits",(start_timeindex+m_MWADataCube.m_StartTimeIndex),dm);
                    if( gDebugLocal >= 1 ){
                        printf("DEBUG1 : szName = %s (m_MWADataCube.m_StartTimeIndex = %d)\n",szName,m_MWADataCube.m_StartTimeIndex);
                    }
 
                    CMWAFits* pDedispersedImage = new CMWAFits( szName, pImage->GetXSize(), pImage->GetYSize() );
                    pDedispersedImage->SetValue( 0 );
                    pDedispersedImage->SetKeys( pImage->GetKeys() );
                    m_DedispersedImages.push_back( pDedispersedImage );
                }
            }
        }else{
           printf("ERROR in code : pImage = NULL in CDedispSearch::Dedisperse (&m_MWADataCube = %p)\n",(&m_MWADataCube));
           exit(-1);
        }
    }
    
    printf("CDedispSearch::Dedisperse(dm=%.2f) : de-dispersing images range %d - %d\n",dm,0,m_MWADataCube.m_Timesteps);
    for(int start_timeindex=0;start_timeindex<m_MWADataCube.m_Timesteps;start_timeindex++){
       CMWAFits* pDedispersedImage = m_DedispersedImages[start_timeindex];
       pDedispersedImage->SetValue( 0 );
    }
    printf("DEBUG : allocted %d buffers for de-dispersed images\n",(int)(m_DedispersedImages.size()));
    m_DedispersedImagesIndex = 0;

    // clear dedispersed images, but do not remove 
    // CleanDedispersed( 0 );
    
    
    double inttime=0.5;
    if( m_MWADataCube.m_Metafits ){
       inttime = m_MWADataCube.m_Metafits->inttime;
    }else{
       printf("WARNING : m_MWADataCube.m_Metafits = NULL -> cannot get real integration time !\n");
    }
    
    // quick and to match python version of the code :
    double freq_upper_end_mhz = m_MWADataCube.m_FreqUpperMHz;    
    double delta_time         = inttime;
    double min_allowed_rms    = 0.00001;
    
    for(int start_timeindex=0;start_timeindex<m_MWADataCube.m_Timesteps;start_timeindex++){
        double unixtime = m_MWADataCube.m_StartUnixTime + start_timeindex * inttime;
        
        CMWAFits* pStartImage = m_MWADataCube.GetImage( 0, start_timeindex );
//        if( !pStartImage ){
//           continue;
//        }
        double start_unixtime = m_MWADataCube.m_StartUnixTime;
        if( pStartImage ){
           start_unixtime = pStartImage->GetUnixTime();
        }
        
        int used_channels = 0;        
        for(int cc=0;cc<m_MWADataCube.m_Channels;cc++){        
            if( cc < m_MWADataCube.m_CoarseChannels.size() ){        
//                double unixtime_fits = GetUnixTime();
                
                int coarse_channel = m_MWADataCube.m_CoarseChannels[cc];                
                double freq_cc_mhz = coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; //  should be added to calculate arrival time at the higher end of the coarse channel
            
                if( gDebugLocal >= 2 ){
                    printf("\tCalculating arrival time from upper end of frequency at %.4f MHz to upper end of channel = %d ( %d / %.4f MHz )\n",freq_upper_end_mhz,cc,coarse_channel,freq_cc_mhz);
                }
            
                double dispersion_delay_ms = dispersion_milisec( freq_upper_end_mhz, freq_cc_mhz, dm ); // # freq_mhz+23.00*1.28+0.5*1.28 - is the higher end of freq. band 
                double dispersion_delay_sec = dispersion_delay_ms / 1000.00;
                double dispersion_delay_timesteps_float = dispersion_delay_sec / delta_time;

// unitl I started to implement multi-obsID version it was int() , but decided to start testing round() :
//                int dispersion_delay_timesteps_int = int( dispersion_delay_timesteps_float );
// Nope I think it was int() on purpose, because we want to include image when the burst started not the next one as in that case it would only 
// containt signal if considerably scattered. 
// SO : it should indeed be int()
// But how to implement int() in m_MWADataCube.GetClosestImage ???
                int dispersion_delay_timesteps_int = int( dispersion_delay_timesteps_float );
            
                // time index :
                int timeindex_estimate = start_timeindex + dispersion_delay_timesteps_int;
                // Proper calculation using time :
                double image_unixtime = start_unixtime + dispersion_delay_sec;
                double min_dist=1e20;
                int timeindex = m_MWADataCube.GetOptimalImage( image_unixtime, min_dist ); // was GetClosestImage
// WORKS OK :                int timeindex = timeindex_estimate;
                if( gDebugLocal >= 2 ){
                   printf("DEBUG : image_unixtime = %.2f ( = %.2f + %.2f ) -> closest timeindex = %d ( estimated guess = %d = %d + %d )\n",image_unixtime,m_MWADataCube.m_StartUnixTime,dispersion_delay_sec,timeindex,timeindex_estimate,start_timeindex , dispersion_delay_timesteps_int);
                }
            
                if( timeindex < m_MWADataCube.m_Timesteps && timeindex >= 0 ){
                    CMWAFits* pImage = m_MWADataCube.GetImage( cc, timeindex );
                    if( pImage ){
                        double image_mean=0.00,image_rms=0.00,image_minval=0.00,image_maxval=0.00;
                        int border = 10;
                        // pImage->GetStatBorder( image_mean, image_rms, image_minval, image_maxval, border );
                        pImage->GetStat( border );
                        image_mean = pImage->m_Mean;
                        image_rms  = pImage->m_RMS;
                        image_minval = pImage->m_Min;
                        image_maxval = pImage->m_Max;
                        if( gDebugLocal >= 2 ){
                            printf("Image( channel=%d , timeindex=%d ) : mean = %.2f Jy, rms = %.2f Jy, min_val = %.2f Jy, max_val = %.2f Jy\n",cc, timeindex, image_mean, image_rms, image_minval, image_maxval );
                        }
                    
                        if ( image_rms > min_allowed_rms && image_rms < m_MaxRMSOnSingle ){
                            // add to dedispersed image:
                            char szName[512];
                            sprintf(szName,"starttimdeindex_%04d_dm%.2f.fits",(start_timeindex+m_MWADataCube.m_StartTimeIndex),dm);
                            if( gDebugLocal >= 1 ){
                                printf("DEBUG2 : szName = %s (m_MWADataCube.m_StartTimeIndex = %d)\n",szName,m_MWADataCube.m_StartTimeIndex);
                            }

                            /*if( m_DedispersedImagesIndex >= m_DedispersedImages.size() ){
                                CMWAFits* pDedispersedImage = new CMWAFits( szName, pImage->GetXSize(), pImage->GetYSize() );                                                                                                                                
                                pDedispersedImage->SetValue( 0 );
                                m_DedispersedImages.push_back( pDedispersedImage );
                                m_DedispersedImagesIndex = m_DedispersedImages.size() - 1; // index set to last image 
                            }else{
                                (*(m_DedispersedImages[m_DedispersedImagesIndex])).SetFileName( szName );
                            }*/
                            
                            (*(m_DedispersedImages[m_DedispersedImagesIndex])).SetFileName( szName );
                            (*(m_DedispersedImages[m_DedispersedImagesIndex])) += (*pImage);
                            if( CDedispSearch::m_bDebugDedispersion || 1 ){
                               (*(m_DedispersedImages[m_DedispersedImagesIndex])).m_OriginalImages.push_back( pImage->GetFileName() );
                            }
                            used_channels = used_channels + 1;
                     
                            // add to mean image :
                            // (*m_MWADataCube.m_pMeanImage) += (*pImage);
                            // (m_MWADataCube.m_pMeanImage)->inc_lines_counter();
                    
                            // flag good image :
                            // m_MWADataCube.m_ImageFlags[cc][timeindex] = 1;
                        }else{
                           if( gDebugLocal>=1){
                              printf("\t : Image( channel=%d , timeindex=%d ) skipped due to rms = %.4f outside range (%.4f , %.4f) !!!\n",cc, timeindex, image_rms, min_allowed_rms, m_MaxRMSOnSingle );
                           }
                        }
                    }else{
                        if( gDebugLocal>=1){
                            printf("\tWARNING : Image( channel=%d , timeindex=%d ) = NULL !!!\n",cc, timeindex );
                        }
                    }                
                }else{                    
                    if( gDebugLocal>=2){
                       printf("\tSKIPPED : cc = %.2f MHz : delay = %.2f [sec] -> delta_timesteps=%.2f=%d >= %d timesteps , timeindex = %d",freq_cc_mhz,dispersion_delay_sec,dispersion_delay_timesteps_float,dispersion_delay_timesteps_int,m_MWADataCube.m_Timesteps,timeindex);                       
                    }
                    if( timeindex < 0 ){
                       printf("\t\tWARNING : timeindex = %d !!!\n",timeindex);
                    }
                }
            }else{
               if( gDebugLocal>=2){printf("WARNING : missing image at coarse channel = %d\n",cc);}
            }
        }        

        // check if number of used channels is larger than required minimum :        
        // if( used_channels > m_MinUsedChannels ){
        // 2019-01-03 : commented in order to save all images (even bad ones, but do not analyse bad ones)
        CMWAFits* pDedispersedImage = (m_DedispersedImages[m_DedispersedImagesIndex]);
        if( used_channels > 0 ){            
            pDedispersedImage->Divide( used_channels );
            pDedispersedImage->GetStat( 20 );            
        }
        pDedispersedImage->m_UsedChannels = used_channels;

        char out_full_path[1024];
        const char* szName = pDedispersedImage->GetFileName();
        sprintf(out_full_path,"%s/dm_%06.2f/%s",CDedispSearch::m_ResultsPath.c_str(),dm,szName);
//        if( !MyFile::DoesFileExist(out_full_path) || CDedispSearch::m_bForceOverwrite ){
        if( TRUE ){ // always overwrite new files in order that sliding window creates files using more averaged files ! This was a bug causing high RMS !!!
            if( (pDedispersedImage->m_OriginalImages).size() > 0 ){
                printf("INFO : saving file %s build from the following %d images (de-dispersion) :\n",out_full_path,int((pDedispersedImage->m_OriginalImages).size()));
            }else{
                printf("INFO : saving file %s\n",out_full_path);
            }
            if( pDedispersedImage->WriteFits( out_full_path ) ){
                printf("ERROR : could not write de-dispersed FITS file %s\n",out_full_path);                                    
            }
            if( gDebugLocal>=2){
                printf("DEBUG : saved file %s with value %.2f at (0,0)\n",out_full_path,pDedispersedImage->getXY(0,0));
            }
            
            if( (pDedispersedImage->m_OriginalImages).size() > 0 ){
               for(int m=0;m<(pDedispersedImage->m_OriginalImages).size();m++){
                  printf("\t%s\n",(pDedispersedImage->m_OriginalImages)[m].c_str());
               }
            }
        }else{
            printf("INFO : de-dispersed file %s already exists\n",out_full_path);
        }        

        m_DedispersedImagesIndex++; // set index to next image (it should really be called startime_index)
        // }
    }
    
}

// here I go over timesteps (not coarse channels to include much more pixels !)
void CDedispSearch::Dedisperse( double dm )
{
    // prepare memory for de-dispersed images :
    if( m_DedispersedImages.size() < m_MWADataCube.m_Timesteps ){
        printf("DEBUG : getting image from timestep = %d , channel = %d\n",m_MWADataCube.m_FirstCorrectTimestep, m_MWADataCube.m_FirstCorrectChannel );
        CMWAFits* pImage = m_MWADataCube.GetImage( m_MWADataCube.m_FirstCorrectChannel , m_MWADataCube.m_FirstCorrectTimestep );
        
        if( pImage ){
            for(int start_timeindex=0;start_timeindex<m_MWADataCube.m_Timesteps;start_timeindex++){
                if( start_timeindex >= m_DedispersedImages.size() ){
                    char szName[512];
                    sprintf(szName,"starttimdeindex_%04d_dm%.2f.fits",(start_timeindex+m_MWADataCube.m_StartTimeIndex),dm);
                    if( gDebugLocal >= 1 ){
                        printf("DEBUG1 : szName = %s (m_MWADataCube.m_StartTimeIndex = %d)\n",szName,m_MWADataCube.m_StartTimeIndex);
                    }
 
                    CMWAFits* pDedispersedImage = new CMWAFits( szName, pImage->GetXSize(), pImage->GetYSize() );
                    pDedispersedImage->SetValue( 0 );
                    pDedispersedImage->SetKeys( pImage->GetKeys() );
                    m_DedispersedImages.push_back( pDedispersedImage );
                }
            }
        }else{
           printf("ERROR in code : pImage = NULL in CDedispSearch::Dedisperse (&m_MWADataCube = %p)\n",(&m_MWADataCube));
           exit(-1);
        }
    }
    
    printf("CDedispSearch::Dedisperse(dm=%.2f) : de-dispersing images range %d - %d\n",dm,0,m_MWADataCube.m_Timesteps);
    for(int start_timeindex=0;start_timeindex<m_MWADataCube.m_Timesteps;start_timeindex++){
       CMWAFits* pDedispersedImage = m_DedispersedImages[start_timeindex];
       pDedispersedImage->SetValue( 0 );
    }
    printf("DEBUG : allocted %d buffers for de-dispersed images\n",(int)(m_DedispersedImages.size()));
    m_DedispersedImagesIndex = 0;

    // clear dedispersed images, but do not remove 
    // CleanDedispersed( 0 );
    
    
    double inttime=0.5;
    if( m_MWADataCube.m_Metafits ){
       inttime = m_MWADataCube.m_Metafits->inttime;
    }else{
       printf("WARNING : m_MWADataCube.m_Metafits = NULL -> cannot get real integration time !\n");
    }
    
    // quick and to match python version of the code :
    double freq_upper_end_mhz = m_MWADataCube.m_FreqUpperMHz;    
    double freq_lower_end_mhz = MWA_COARSE_CHANNEL/2.00 + m_MWADataCube.m_CoarseChannels[0]*MWA_COARSE_CHANNEL;
    double delta_time         = inttime;
    double min_allowed_rms    = 0.00001;
    
    for(int start_timeindex=0;start_timeindex<m_MWADataCube.m_Timesteps;start_timeindex++){
        double unixtime = m_MWADataCube.m_StartUnixTime + start_timeindex * inttime;
        
        CMWAFits* pStartImage = m_MWADataCube.GetImage( 0, start_timeindex );
//        if( !pStartImage ){
//           continue;
//        }
        double start_unixtime = m_MWADataCube.m_StartUnixTime;
        if( pStartImage ){
           start_unixtime = pStartImage->GetUnixTime();
        }
        
        int used_channels = 0;        
        double arrival_time = 0.00 , arrival_time_step = 0.50;
        double freq_cc_mhz = freq_upper_end_mhz;
        int timeindex = -1;
        double min_dist = 0.00;
        
        // this loop goes over arrival times 
        while( freq_cc_mhz >= freq_lower_end_mhz && arrival_time<=(inttime*m_MWADataCube.m_Timesteps) && timeindex < m_MWADataCube.m_Timesteps ){
            double freq_mhz = freq_cc_mhz;
            double coarse_channel = freq_mhz / MWA_COARSE_CHANNEL; // floor -> absolute number of coarse channel 
            
            // find index in the list of coarse channels in the metadata :
            int cc = find_value_double( m_MWADataCube.m_CoarseChannels, coarse_channel ); // coarse_channel - m_MWADataCube.m_CoarseChannels[0]; // WARNING : assumes it's not a picket fence observation !
            
            if( cc < 0 ){
               printf("ERROR : freq_mhz = %.2f MHz -> coarse channel = %.6f not found in the observation metadata - this means bug in the code (as it should not happen) -> cannot continue / exiting function\n",freq_mhz,coarse_channel);
               return;                
            }
            
            if( cc < m_MWADataCube.m_CoarseChannels.size() ){        
//                double unixtime_fits = GetUnixTime();                
                if( gDebugLocal >= 2 ){
                    printf("\tCalculating arrival time from upper end of frequency at %.4f MHz to upper end of channel = %d ( %.6f / %.4f MHz )\n",freq_upper_end_mhz,cc,coarse_channel,freq_cc_mhz);
                }
            
                double dispersion_delay_ms = arrival_time*1000; //  dispersion_milisec( freq_upper_end_mhz, freq_cc_mhz, dm ); // # freq_mhz+23.00*1.28+0.5*1.28 - is the higher end of freq. band 
                double dispersion_delay_sec = dispersion_delay_ms / 1000.00;
                double dispersion_delay_timesteps_float = dispersion_delay_sec / delta_time;

// unitl I started to implement multi-obsID version it was int() , but decided to start testing round() :
//                int dispersion_delay_timesteps_int = int( dispersion_delay_timesteps_float );
// Nope I think it was int() on purpose, because we want to include image when the burst started not the next one as in that case it would only 
// containt signal if considerably scattered. 
// SO : it should indeed be int()
// But how to implement int() in m_MWADataCube.GetClosestImage ???
                int dispersion_delay_timesteps_int = int( dispersion_delay_timesteps_float );
            
                // time index :
                int timeindex_estimate = start_timeindex + dispersion_delay_timesteps_int;
                // Proper calculation using time :
                double image_unixtime = start_unixtime + dispersion_delay_sec;
                timeindex = m_MWADataCube.GetOptimalImage( image_unixtime , min_dist ); // was GetClosestImage
// WORKS OK :                int timeindex = timeindex_estimate;
                if( gDebugLocal >= 2 ){
                   printf("DEBUG : image_unixtime = %.2f ( = %.2f + %.2f ) -> closest timeindex = %d ( estimated guess = %d = %d + %d ) , min_dist = %.2f [sec]\n",image_unixtime,m_MWADataCube.m_StartUnixTime,dispersion_delay_sec,timeindex,timeindex_estimate,start_timeindex , dispersion_delay_timesteps_int,min_dist);
                }
            
                if( timeindex < m_MWADataCube.m_Timesteps && timeindex >= 0 && min_dist <= inttime ){
                    CMWAFits* pImage = m_MWADataCube.GetImage( cc, timeindex );
                    if( pImage ){
                        double image_mean=0.00,image_rms=0.00,image_minval=0.00,image_maxval=0.00;
                        int border = 10;
                        // pImage->GetStatBorder( image_mean, image_rms, image_minval, image_maxval, border );
                        pImage->GetStat( border );
                        image_mean = pImage->m_Mean;
                        image_rms  = pImage->m_RMS;
                        image_minval = pImage->m_Min;
                        image_maxval = pImage->m_Max;
                        if( gDebugLocal >= 2 ){
                            printf("Image( channel=%d , timeindex=%d ) : mean = %.2f Jy, rms = %.2f Jy, min_val = %.2f Jy, max_val = %.2f Jy\n",cc, timeindex, image_mean, image_rms, image_minval, image_maxval );
                        }
                    
                        if ( image_rms > min_allowed_rms && image_rms < m_MaxRMSOnSingle ){
                            // add to dedispersed image:
                            char szName[512];
                            sprintf(szName,"starttimdeindex_%04d_dm%.2f.fits",(start_timeindex+m_MWADataCube.m_StartTimeIndex),dm);
                            if( gDebugLocal >= 1 ){
                                printf("DEBUG2 : szName = %s (m_MWADataCube.m_StartTimeIndex = %d)\n",szName,m_MWADataCube.m_StartTimeIndex);
                            }

                            /*if( m_DedispersedImagesIndex >= m_DedispersedImages.size() ){
                                CMWAFits* pDedispersedImage = new CMWAFits( szName, pImage->GetXSize(), pImage->GetYSize() );                                                                                                                                
                                pDedispersedImage->SetValue( 0 );
                                m_DedispersedImages.push_back( pDedispersedImage );
                                m_DedispersedImagesIndex = m_DedispersedImages.size() - 1; // index set to last image 
                            }else{
                                (*(m_DedispersedImages[m_DedispersedImagesIndex])).SetFileName( szName );
                            }*/
                            
                            (*(m_DedispersedImages[m_DedispersedImagesIndex])).SetFileName( szName );
                            (*(m_DedispersedImages[m_DedispersedImagesIndex])) += (*pImage);
                            if( CDedispSearch::m_bDebugDedispersion || 1 ){
                               (*(m_DedispersedImages[m_DedispersedImagesIndex])).m_OriginalImages.push_back( pImage->GetFileName() );
                            }
                            used_channels = used_channels + 1;
                     
                            // add to mean image :
                            // (*m_MWADataCube.m_pMeanImage) += (*pImage);
                            // (m_MWADataCube.m_pMeanImage)->inc_lines_counter();
                    
                            // flag good image :
                            // m_MWADataCube.m_ImageFlags[cc][timeindex] = 1;
                        }else{
                           if( gDebugLocal>=1){
                              printf("\t : Image( channel=%d , timeindex=%d ) skipped due to rms = %.4f outside range (%.4f , %.4f) !!!\n",cc, timeindex, image_rms, min_allowed_rms, m_MaxRMSOnSingle );
                           }
                        }
                    }else{
                        if( gDebugLocal>=1){
                            printf("\tWARNING : Image( channel=%d , timeindex=%d ) = NULL !!!\n",cc, timeindex );
                        }
                    }                
                }else{                    
                    if( gDebugLocal>=2){
                       printf("\tSKIPPED : cc = %.2f MHz : delay = %.2f [sec] -> delta_timesteps=%.2f=%d >= %d timesteps , timeindex = %d , min_dist = %.4f [sec]",freq_cc_mhz,dispersion_delay_sec,dispersion_delay_timesteps_float,dispersion_delay_timesteps_int,m_MWADataCube.m_Timesteps,timeindex,min_dist);                       
                    }
                    if( timeindex < 0 ){
                       printf("\t\tWARNING : timeindex = %d !!!\n",timeindex);
                    }
                }
            }else{
               if( gDebugLocal>=2){printf("WARNING : missing image at coarse channel = %d\n",cc);}
            }
            
            arrival_time += arrival_time_step;
            freq_cc_mhz = dispersion_freq_mhz( freq_upper_end_mhz , arrival_time*1000.00, dm );
        }        

        // check if number of used channels is larger than required minimum :        
        // if( used_channels > m_MinUsedChannels ){
        // 2019-01-03 : commented in order to save all images (even bad ones, but do not analyse bad ones)
        CMWAFits* pDedispersedImage = (m_DedispersedImages[m_DedispersedImagesIndex]);
        if( used_channels > 0 ){            
            pDedispersedImage->Divide( used_channels );
            pDedispersedImage->GetStat( 20 );            
        }
        pDedispersedImage->m_UsedChannels = used_channels;

        char out_full_path[1024];
        const char* szName = pDedispersedImage->GetFileName();
        sprintf(out_full_path,"%s/dm_%06.2f/%s",CDedispSearch::m_ResultsPath.c_str(),dm,szName);
//        if( !MyFile::DoesFileExist(out_full_path) || CDedispSearch::m_bForceOverwrite ){
        if( TRUE ){ // always overwrite new files in order that sliding window creates files using more averaged files ! This was a bug causing high RMS !!!
            if( (pDedispersedImage->m_OriginalImages).size() > 0 ){
                printf("INFO : saving file %s build from the following %d images (de-dispersion) :\n",out_full_path,int((pDedispersedImage->m_OriginalImages).size()));
            }else{
                printf("INFO : saving file %s\n",out_full_path);
            }
            if( pDedispersedImage->WriteFits( out_full_path ) ){
                printf("ERROR : could not write de-dispersed FITS file %s\n",out_full_path);                                    
            }
            if( gDebugLocal>=2){
                printf("DEBUG : saved file %s with value %.2f at (0,0)\n",out_full_path,pDedispersedImage->getXY(0,0));
            }
            
            if( (pDedispersedImage->m_OriginalImages).size() > 0 ){
               for(int m=0;m<(pDedispersedImage->m_OriginalImages).size();m++){
                  printf("\t%s\n",(pDedispersedImage->m_OriginalImages)[m].c_str());
               }
            }
        }else{
            printf("INFO : de-dispersed file %s already exists\n",out_full_path);
        }        

        m_DedispersedImagesIndex++; // set index to next image (it should really be called startime_index)
        // }
    }
    
}

double CDedispSearch::InjectFRB( CMWAFits& dynspec, int start_timeindex, double value, bool overwrite /* =false */, double snr /*=-100*/, int* end_timeindex /*=NULL*/ )
{
    if( !m_MWADataCube.m_Metafits ){
       printf("ERROR : to correctly generate time series images metafits file is required !\n");
       exit(-1);
    }

    if( end_timeindex ){
       (*end_timeindex) = start_timeindex;
    }
 
    double return_value = 1.00;
    double inttime = m_MWADataCube.m_Metafits->inttime;
    

    int n_coarse_channels = m_MWADataCube.m_CoarseChannels.size();
    int first_coarse_channel  = m_MWADataCube.m_CoarseChannels[0];
    double freq_lower_end_mhz = first_coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // MWA_COARSE_CHANNEL/2.00 + m_MWADataCube.m_CoarseChannels[0]*MWA_COARSE_CHANNEL;
    double freq_upper_end_mhz = m_MWADataCube.m_CoarseChannels[n_coarse_channels-1]*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // m_MWADataCube.m_FreqUpperMHz;
    double delta_time = 0.1; // 0.001 too much memory ~14 GB , 0.01 - 1.4 GB but too long X-axis !

    int dm_size = ( (m_MaxDM - m_MinDM) / m_StepDM ) + 1;
   
    m_MaxDispersiveSweep = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, m_MaxDM )/1000.00; // in seconds
    int n_time_steps = dynspec.GetXSize();
    
    // go trough all DM templetas and add pixels in these templates to calculate sum.
//    sort( m_DispersionTemplates.begin(), m_DispersionTemplates.end() );
    for(map< string, CDispersedPixels >::iterator it = m_DispersionTemplates.begin(); it != m_DispersionTemplates.end(); ++it) {
//    for(int ii=0;ii<m_DispersionTemplates.size(); ii++){
//         CDispersedPixels>::iterator* it = &( m_DispersionTemplates[ii] );
      
         CDispersedPixels& pixels = it->second;
         double dm = atof( it->first.c_str() );
         double sweep_time_dm = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, dm )/1000.00; // in seconds
      
         double start_uxtime = -sweep_time_dm;
         double end_uxtime   = dynspec.GetXSize()*inttime + sweep_time_dm;
         printf("Adding FRB with dm = %.2f\n",dm);
         
         // snr parameter overwrites value:         
         double mean_local = 0.00, rms_local = 0.00;
         if( snr > 0 ){
            // calculate local RMS and derive value as SNR*RMS :
            int n_timesteps_for_rms=50;
            double sum=0;
            double sum2=0;
            int    sum_cnt=0;
            for(int t=(start_timeindex-n_timesteps_for_rms);t<=(start_timeindex+n_timesteps_for_rms);t++){
              if( t>=0 && t<dynspec.GetXSize() ){
                 double sum_signal = 0.00;
                 bool   ok = true;
                 for(int p=0;p<pixels.size();p++){
                    CPixel& pixel = pixels[p];
                    int timeindex = t + pixel.time;
                    int cc = pixel.cc;
                    
                    if( timeindex <= dynspec.GetXSize() ){
                       sum_signal += dynspec.getXY( timeindex , cc );
                    }else{
                       ok = false;
                    }
                 }
                 if( ok ){
                    sum += sum_signal;
                    sum2 += sum_signal*sum_signal;
                    sum_cnt++;
                 }
              } 
            }            
            
            mean_local = sum/sum_cnt;
            rms_local = sqrt( sum2/sum_cnt - mean_local*mean_local );
            value = snr*rms_local / pixels.size();
            return_value = snr*rms_local;            
            printf("DEBUG : snr=%.4f and mean and rms local = %.4f and %.4f -> value := %.4f same for all %d pixels ( total fluence = %.4f )\n",snr,mean_local,rms_local,value,int(pixels.size()),(snr*rms_local));
         }

         double sum_check = 0.00;         
         for(int p=0;p<pixels.size();p++){
               CPixel& pixel = pixels[p];
               int timeindex = start_timeindex + pixel.time;
               int cc = pixel.cc;
               
               if( end_timeindex ){
                  if( timeindex > (*end_timeindex) ){
                     (*end_timeindex) = timeindex;
                  }
               }

               if( overwrite ){               
                  // overwrite whatever is there :
                  dynspec.setXY( timeindex, cc, value*pixel.weight );
                  sum_check += value;
               }else{
                  // more realistic add to existing noise :
                  double current_value = dynspec.getXY( timeindex, cc );
                  double new_value =  current_value + value*pixel.weight;
                  printf("(%d,%d) : %.4f + %.4f = %.4f Jy\n", timeindex, cc, current_value, value, new_value);
                  dynspec.setXY( timeindex, cc, new_value );
                  
                  sum_check += new_value;
               }
         }         
         
         if( snr > 0 ){
            printf("DEBUG : injection of SNR=%.4f signal ended up with total_sum = %.4f which is snr = %.4f\n",snr,sum_check,(sum_check-mean_local)/rms_local);
         }
         
         break;
     }


     return return_value;
}

// rms_map_new is map of RMS-s calculated as Sum of RMS_channel^2 for all pixels in the DM-sweep 
int CDedispSearch::GetDedispersedSeries( CMWAFits& dynspec, CMWAFits& dedispersed_series, CMWAFits& out_count_map, CMWAFits& rms_map_new,
                                         std::vector<double>* p_channels_list, double threshold_in_sigma, 
                                         MyOFile* candidates_file, int max_radius, double LogBelowDM /*=-1e6*/ , bool bIgnoreMaxValue /* = false */,
                                         bool bDoWeigthByRMS  )
{

    // temporary :
/*    for(int y=0;y<dynspec.GetYSize();y++){
       for(int x=0;x<dynspec.GetXSize();x++){
          dynspec.setXY( x, y, CMyFit::GetGaussFast( 2.773, 0.0296 ) );
       }
    }*/

// bDoWeigthByRMS = false;   

    if( !m_MWADataCube.m_Metafits ){
       printf("ERROR : to correctly generate time series images metafits file is required !\n");
       exit(-1);
    }

    // calculate mean and rms in a row to weight pixels :
    CBgArray mean_channel,rms_channel;
    dynspec.MeanLines( mean_channel , rms_channel );
    // save to text files named accoring to name of input dynamic spectrum:
    string szOutTmpFile;
    change_ext( dynspec.GetFileName(), "_rms_channel.txt", szOutTmpFile, true );
    rms_channel.SaveToFile( szOutTmpFile.c_str() );
    change_ext( dynspec.GetFileName(), "_mean_channel.txt", szOutTmpFile, true );
    mean_channel.SaveToFile( szOutTmpFile.c_str() );
    
    double inttime = m_MWADataCube.m_Metafits->inttime;
    

    int n_coarse_channels = m_MWADataCube.m_CoarseChannels.size();
    int first_coarse_channel  = m_MWADataCube.m_CoarseChannels[0];
    double freq_lower_end_mhz = first_coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // MWA_COARSE_CHANNEL/2.00 + m_MWADataCube.m_CoarseChannels[0]*MWA_COARSE_CHANNEL;
    double freq_upper_end_mhz = m_MWADataCube.m_CoarseChannels[n_coarse_channels-1]*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // m_MWADataCube.m_FreqUpperMHz;
    double delta_time = 0.1; // 0.001 too much memory ~14 GB , 0.01 - 1.4 GB but too long X-axis !
    if( delta_time > inttime ){
       delta_time = inttime; // /10.00;
    }
    double delta_time_step = delta_time / 4.00; // this is just a step to go through the time

    int dm_size = ( (m_MaxDM - m_MinDM) / m_StepDM ); // + 1; 2023-07-17 removed +1 
   
    m_MaxDispersiveSweep = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, m_MaxDM )/1000.00; // MAX DM sweep time in seconds
    // m_MaxDispersiveSweep = max_sweep_time_dm;
    double max_sweep_time_dm = 0.00; // this is set to non-zero value when m_bAllowPreStartArrivals=true -> meaning that we want to calculate DTS and include pre-start data arrival times
    double total_duration = dynspec.GetXSize()*inttime;
    if( m_bAllowPreStartArrivals ){
       max_sweep_time_dm = m_MaxDispersiveSweep;
       total_duration = (max_sweep_time_dm + dynspec.GetXSize()*inttime); // was + max_sweep_time_dm); - but this is wrong as any FRB which arrived at the upper end of freq. band after the end of          
                                                           // the observation will not be detected - Even a single pixel
    }                                                           
    int timeseries_size = ( total_duration / delta_time ) + 1;
    int n_time_steps = dynspec.GetXSize();
    
    dedispersed_series.Realloc( timeseries_size, dm_size );
    dedispersed_series.SetValue(0.00);
    
    out_count_map.Realloc( timeseries_size, dm_size );
    out_count_map.SetValue(0.00);
    
    rms_map_new.Realloc( timeseries_size, dm_size );
    rms_map_new.SetValue(0.00);
    
    printf("DEBUG : DTS size %d x %d , dynaspec size %d x %d\n",timeseries_size,dm_size,dynspec.GetXSize(),dynspec.GetYSize());
    
    // ONE TIME TEST TO SAVE the unique paths. A unique path is DM+StartTimeIndex, now I have a code below which counts them always - no need to change 
    // the code for every test !
/*    static bool bCountedPaths=true;  // do not write file
    FILE* out_f = NULL;
    if( !bCountedPaths ){
       out_f = fopen("sweep_paths.txt","w");
       bCountedPaths = true;
    }*/
    int max_time_index = dynspec.GetXSize()*2;
    int n_unique_sweep_paths = 0;
    int template_index = 0;
    vector<double> dm_values;
    int max_count=-1; // return value - maximum number of pixels added 

    // go trough all DM templetas and add pixels in these templates to calculate sum.
//    sort( m_DispersionTemplates.begin(), m_DispersionTemplates.end() );
    for(map< string, CDispersedPixels >::iterator it = m_DispersionTemplates.begin(); it != m_DispersionTemplates.end(); ++it) {
//    for(int ii=0;ii<m_DispersionTemplates.size(); ii++){
//         CDispersedPixels>::iterator* it = &( m_DispersionTemplates[ii] );
      
      
         vector<int> path_counter; // count starttimes for given DM
         path_counter.assign( max_time_index, 0 );
         
         CDispersedPixels& pixels = it->second;
         double dm = atof( it->first.c_str() );
         double sweep_time_dm = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, dm )/1000.00; // in seconds

         double start_uxtime = 0.00;
         if( m_bAllowPreStartArrivals ){      
            start_uxtime = -sweep_time_dm;  // start unix time for a given dm (not max DM)
         }
         double end_uxtime   = dynspec.GetXSize()*inttime; // was  + sweep_time_dm; - but this does not make sense as any FRB which arrived at the upper end of freq. band after the end of 
                                                           // the observation will not be detected - Even a single pixel
         printf("\nAnalysing dm = %.2f (%d pixels in dispersion path , template_index = %d)\n",dm,int(pixels.size()),template_index);
         
         while( start_uxtime <= end_uxtime ){
            // 2023-07-15 : this "-1" was causing first column to be the same in DTS : see page 3 in /home/msok/Desktop/PAWSEY/PaCER/logbook/20230715_changes_in_frbsearch_code.odt
            //     this is a big change so I've keep the previous version in : /home/msok/bighorns/software/analysis/frb_search/backup/20230715/1115am/ and /home/msok/github/msfitslib/apps/frb_search/backup/20230715/1117am
            // int start_timeindex = (start_uxtime / inttime) - 1;            
            // 2023-07-15 : it should really be uncommented but I will minimise changes for now:
            int start_timeindex = int(start_uxtime / inttime);
            //if( start_timeindex < 0 ){
            //   start_timeindex = 0; // 2023-07-15 : do not allow negative values of this index !!! I think this was doubling same values at time=0 and time=1 in DTS
            //}
            if( start_timeindex < path_counter.size() ){
               if( dm >= 200 ){ // only count DM >= 200 
                  path_counter[start_timeindex] = path_counter[start_timeindex] + 1;
               }
            }
            
/* DEBUG            
            if( (start_timeindex >= 3428 && dm>500) || start_timeindex==-8 ){
               printf("odo\n");
            } */
         
            double sum = 0 , sum_weight = 0.00 , max_value = -1e20, max_weight = -1e20, sum_sigma2 = 0.00;
            int    count = 0 , count_with_nan = 0;
            if( dm < LogBelowDM ){            
               printf("\tDUMP_PIxELS : DM = %.2f, start_uxtime = %.2f : ",dm,start_uxtime);
            }
//            vector<double> dm_values;
            if( dm_values.size() > 0 ){
               dm_values.clear();
            }
            for(int p=0;p<pixels.size();p++){
               CPixel& pixel = pixels[p];
               int timeindex = start_timeindex + pixel.time;
               
               // TODO : add if pixel.cc in range of RFI flagged - add member variable or parameter here !
               bool is_chan_flagged = dynspec.is_channel_flagged( pixel.cc );
               if( is_chan_flagged ){
                  continue;
               }
            
               if( timeindex >=0 && timeindex <= n_time_steps ){
                   if( pixel.cc >=0 && pixel.cc < n_coarse_channels ){                    
                      double val = dynspec.getXY( timeindex, pixel.cc );
                      
//                      if( true ){
                         // test of independent values 
//                         val = CMyFit::GetGaussFast( 2.773, 0.0296 );
//                         printf("GAUSS = %.8f\n",val);
//                      }
                      
                      if( isnan(val) == 0 ){
                         // skip NaNs !
                         
                         // new weigthing :
                         double rms_ch = rms_channel[pixel.cc];
                         double weight = 1.00;
                         if( bDoWeigthByRMS ){
                            weight = (1.00/rms_ch);
                         }
                                                  
                         double new_val = val*weight;                         
                         sum += new_val;                  
                         sum_weight += weight;
                         sum_sigma2 += rms_ch*rms_ch; // sum of RMS^2 (calculated for a given channel) , https://en.wikipedia.org/wiki/Sum_of_normally_distributed_random_variables
                         
                         if( m_DebugDM > 0 && fabs(m_DebugDM-dm)<0.001 ){
                            dm_values.push_back( new_val );
                         }
                         
                         if( (val*weight) > max_value ){
                            max_value = (val*weight);
                            max_weight = weight;
                         }
                         
                         if( dm < LogBelowDM ){
                            printf("(%d,%d)=%.4f->%.4f ,",timeindex,pixel.cc,val,val*weight);
                         }
                         
//                         if( out_f ){
//                            all_sweep_paths.push_back( pixels );
//                            fprintf(out_f,"%s-%08d\n",it->first.c_str(),start_timeindex);
//                         }
                              
                         count++;     
                      }
                      count_with_nan++;
                   }else{
                      if( CDedispSearch::m_bAllowNegativeDM ){
                         ; // all good when allowing for negative DMs 
                      }else{
                         printf("ERROR in code GetDedispersedSeries - cc <0 or c>=%d !!! -> cannot continue !!!\n",n_coarse_channels);
                         exit(-1);                   
                      }
                   }
               }
            }
            double sum_weighted = sum / sum_weight; // normally sum is divded by sum_weight which when RMS-weigiting option (-w 0) is not used is sum_weight=COUNT (so it is already per pixel)
            if( count <= 0 && isnan(sum_weighted) ){
               sum_weighted = 0.00;
            }
            
            if( bIgnoreMaxValue ){
               // in this case substract maximum value (to avoid 2 loops) and its weight (Note : max_weight is weight corresponding to max_value not MAX(weights) ! )
               sum_weighted = ( sum - max_value ) / ( sum_weight - max_weight );
               count--; // decrease counter
            }
            
            // 2023-07-15 : round(0.5) = 1 -> repetition of the first column in DTS !!! - checking int() -> did not really help still first column is the same ... - leaving it for now 
            int out_timeindex = round((start_uxtime - (-fabs(max_sweep_time_dm)))/delta_time); // this is index in image FITS file (must start from 0 corresponding to -T(DM_max) ,  WAS WRONG (start_uxtime + sweep_time_dm)/delta_time;
                                // int( ( sweep_time_dm - max_sweep_time_dm ) / inttime );
            // 2023-02-22 : was int() but this causes some DM rows to be ZERO !!!                                
//            int out_dm_index = round( (dm-m_MinDM) / m_StepDM); // 2023-02-22 : was int() but this causes some DM rows to be ZERO !!!
            // 2023-07-17 : when changed in CalcTemplatesExact have to also change here:
            int out_dm_index = round( (dm-m_MinDM-m_StepDM/2.00) / m_StepDM);

            if( dm < LogBelowDM ){
               printf("\t -> %.4f (start at pixel (%d,%d) )\n",(sum_weighted),out_timeindex,out_dm_index);fflush(stdout);
            }

            
//            double sum_weighted = sum * ( pixels.size() / count );
            if( count_with_nan <= 0 ){
               printf("\tWARNING : count(including nan)=0 (count non-NaN = %d) for dm=%.2f (sweep time = %.6f [sec]), start_uxtime=%.6f (start_timeindex = %d ) vs. end_uxtime = %.6f\n",count,dm,sweep_time_dm,start_uxtime,start_timeindex,end_uxtime);
            }

//            double sum_weighted = sum / count; // warning first start times with start_uxtime<0 give count=0 -> NaN , but I leave it like this for now ...
            dedispersed_series.setXY( out_timeindex, out_dm_index, sum_weighted );
            out_count_map.setXY( out_timeindex, out_dm_index, count );
            
            double rms_new = sqrt( sum_sigma2 )/count;
            rms_map_new.setXY( out_timeindex, out_dm_index, rms_new );
            if( count > max_count ){
               max_count = count;
            }
            
            if( dm_values.size() ){
               printf("DEBUG : start_ux = %.4f : ",start_uxtime);
               for(int ii=0;ii<dm_values.size();ii++){
                  printf("%.4f + ",dm_values[ii]);
               }
               printf(" = %.6f\n",sum_weighted);
            }
            
            start_uxtime += delta_time_step;
         }
         
         // calculate unique paths:
         for(int pp=0;pp<path_counter.size();pp++){
            if( path_counter[pp] > 0 ){
               n_unique_sweep_paths++;
            }
         }
         
         template_index++;
    }
    
    printf("STAT_INFO : number of unique sweep paths = %d\n",n_unique_sweep_paths);   
    
    dedispersed_series.SetKeyword( "CTYPE2", "DM" );
    dedispersed_series.SetKeyword( "CUNIT2", "1.00" );
    dedispersed_series.SetKeywordFloat( "CRPIX2", 1 );
    dedispersed_series.SetKeywordFloat( "CDELT2", (float)m_StepDM );
    dedispersed_series.SetKeywordFloat( "CRVAL2", (float)(m_MinDM + m_StepDM/2.00) ); // 2023-07-15 : add + m_StepDM/2.00 so that center pixel is 1/2 of DM step (see /home/msok/Desktop/PAWSEY/PaCER/logbook/20230715_changes_in_frbsearch_code.odt )
    dedispersed_series.SetKeyword( "CTYPE1","Time" );
    dedispersed_series.SetKeyword( "CUNIT1","sec");
    dedispersed_series.SetKeyword( "CRPIX1", 1 );
    dedispersed_series.SetKeywordFloat( "CDELT1", delta_time );
    dedispersed_series.SetKeywordFloat( "CRVAL1", -max_sweep_time_dm ); // time starts from 0.00
    dedispersed_series.SetKeyword( "VERSION" , DEDISPERSION_VERSION );
    dedispersed_series.SetKeyword( "USWEEPS", n_unique_sweep_paths );
    dedispersed_series.inttime = delta_time;
    dedispersed_series.m_StartTime = -max_sweep_time_dm;

    out_count_map.SetKeyword( "CTYPE2", "DM" );
    out_count_map.SetKeyword( "CUNIT2", "1.00" );
    out_count_map.SetKeywordFloat( "CRPIX2", 1 );
    out_count_map.SetKeywordFloat( "CDELT2", (float)m_StepDM );
    out_count_map.SetKeywordFloat( "CRVAL2", (float)(m_MinDM + m_StepDM/2.00) ); // 2023-07-15 : add + m_StepDM/2.00 so that center pixel is 1/2 of DM step (see /home/msok/Desktop/PAWSEY/PaCER/logbook/20230715_changes_in_frbsearch_code.odt )
    out_count_map.SetKeyword( "CTYPE1","Time" );
    out_count_map.SetKeyword( "CUNIT1","sec");
    out_count_map.SetKeyword( "CRPIX1", 1 );
    out_count_map.SetKeywordFloat( "CDELT1", delta_time );
    out_count_map.SetKeyword( "CRVAL1", -max_sweep_time_dm ); // time starts from 0.00
    out_count_map.SetKeyword( "VERSION" , DEDISPERSION_VERSION );    
    out_count_map.SetKeyword( "USWEEPS", n_unique_sweep_paths );

    rms_map_new.SetKeyword( "CTYPE2", "DM" );
    rms_map_new.SetKeyword( "CUNIT2", "1.00" );
    rms_map_new.SetKeywordFloat( "CRPIX2", 1 );
    rms_map_new.SetKeywordFloat( "CDELT2", (float)m_StepDM );
    rms_map_new.SetKeywordFloat( "CRVAL2", (float)(m_MinDM + m_StepDM/2.00) ); // 2023-07-15 : add + m_StepDM/2.00 so that center pixel is 1/2 of DM step (see /home/msok/Desktop/PAWSEY/PaCER/logbook/20230715_changes_in_frbsearch_code.odt )
    rms_map_new.SetKeyword( "CTYPE1","Time" );
    rms_map_new.SetKeyword( "CUNIT1","sec");
    rms_map_new.SetKeyword( "CRPIX1", 1 );
    rms_map_new.SetKeywordFloat( "CDELT1", delta_time );
    rms_map_new.SetKeyword( "CRVAL1", -max_sweep_time_dm ); // time starts from 0.00
    rms_map_new.SetKeyword( "VERSION" , DEDISPERSION_VERSION );    
    rms_map_new.SetKeyword( "USWEEPS", n_unique_sweep_paths );

/*    bCountedPaths=true;
    if( out_f ){
       fclose( out_f );
       printf("DEBUG : closed sweep_path.txt file - should not write a new one !\n");
    }*/
    
    return max_count; // maximum number of summed pixels : TODO : do not use too fine delta_time in CalculateTemplatesExact !!!
}


int CDedispSearch::ShowSweepPixels( CMWAFits& dynspec, int start_timeindex,  MyOFile* candidates_file, std::vector<double>* p_channels_list )
{
    if( !m_MWADataCube.m_Metafits ){
       printf("ERROR : to correctly generate time series images metafits file is required !\n");
       exit(-1);
    }

    double inttime = m_MWADataCube.m_Metafits->inttime;
    

    int n_coarse_channels = m_MWADataCube.m_CoarseChannels.size();
    int first_coarse_channel  = m_MWADataCube.m_CoarseChannels[0];
    double freq_lower_end_mhz = first_coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // MWA_COARSE_CHANNEL/2.00 + m_MWADataCube.m_CoarseChannels[0]*MWA_COARSE_CHANNEL;
    double freq_upper_end_mhz = m_MWADataCube.m_CoarseChannels[n_coarse_channels-1]*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // m_MWADataCube.m_FreqUpperMHz;
    double delta_time = 0.1; // 0.001 too much memory ~14 GB , 0.01 - 1.4 GB but too long X-axis !

    int dm_size = ( (m_MaxDM - m_MinDM) / m_StepDM ) + 1;

    double max_sweep_time_dm = 0.00;   
    double total_duration = dynspec.GetXSize()*inttime;  
    if( m_bAllowPreStartArrivals ){
       double max_sweep_time_dm = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, m_MaxDM )/1000.00; // MAX DM sweep time in seconds
       double total_duration = (max_sweep_time_dm + dynspec.GetXSize()*inttime); // was + max_sweep_time_dm); - but this is wrong as any FRB which arrived at the upper end of freq. band after the end of 
                                                                                 // the observation will not be detected - Even a single pixel
    }                                                           
    int timeseries_size = ( total_duration / delta_time ) + 1;
    int n_time_steps = dynspec.GetXSize();
    

    // go trough all DM templetas and add pixels in these templates to calculate sum.
//    sort( m_DispersionTemplates.begin(), m_DispersionTemplates.end() );
    for(map< string, CDispersedPixels >::iterator it = m_DispersionTemplates.begin(); it != m_DispersionTemplates.end(); ++it) {
//    for(int ii=0;ii<m_DispersionTemplates.size(); ii++){
//         CDispersedPixels>::iterator* it = &( m_DispersionTemplates[ii] );
      
         CDispersedPixels& pixels = it->second;
         double dm = atof( it->first.c_str() );
         double sweep_time_dm = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, dm )/1000.00; // in seconds      
         printf("Analysing dm = %.2f\n",dm);
         
/* DEBUG            
            if( (start_timeindex >= 3428 && dm>500) || start_timeindex==-8 ){
               printf("odo\n");
            } */
       
          double sum = 0;
          int    count = 0 , count_with_nan = 0;
          for(int p=0;p<pixels.size();p++){
             CPixel& pixel = pixels[p];
             int timeindex = start_timeindex + pixel.time;
            
             if( timeindex >=0 && timeindex <= n_time_steps ){
                 if( pixel.cc >=0 && pixel.cc < n_coarse_channels ){                    
                    double val = dynspec.getXY( timeindex, pixel.cc );
                    if( isnan(val) == 0 ){
                       // skip NaNs !
                       sum += val;       
                       count++;     
                    }
                    printf("\t(%d,%d) = %.4f -> sum = %.4f\n",timeindex,pixel.cc,val,sum);
                    count_with_nan++;
                 }else{
                    if( CDedispSearch::m_bAllowNegativeDM ){ 
                       ; // all good when allowing for negative DMs 
                    }else{
                       printf("ERROR in code GetDedispersedSeries - cc <0 or c>=%d !!! -> cannot continue !!!\n",n_coarse_channels);
                       exit(-1);                   
                    }
                 }
             }
          }

//          int out_timeindex = int((start_uxtime - (-fabs(max_sweep_time_dm)))/delta_time); // this is index in image FITS file (must start from 0 corresponding to -T(DM_max) ,  WAS WRONG (start_uxtime + sweep_time_dm)/delta_time;
                                // int( ( sweep_time_dm - max_sweep_time_dm ) / inttime );
//          int out_dm_index = int( (dm-m_MinDM) / m_StepDM);
            
//            double sum_weighted = sum * ( pixels.size() / count );
          if( count_with_nan <= 0 ){
             printf("WARNING : count(including nan)=0\n");
          }

          double sum_weighted = sum / count; // warning first start times with start_uxtime<0 give count=0 -> NaN , but I leave it like this for now ...
          printf("DM = %.2f, start_timeindex = %d -> %d pixels give mean = %.4f\n",dm,start_timeindex,count,sum_weighted);
          
          return sum_weighted;
    }
    

    return 1;
}



int CDedispSearch::AnalyseDynSpec( CMWAFits& dynspec, int x, int y, double threshold_in_sigma, MyOFile* candidates_file, int max_radius,  int min_good_pixels, std::vector<double>* p_channels_list  )
{
    char out_full_path[1024];
    sprintf(out_full_path,"%s/dynamic_spectrum_search/%04d_%04d.fits",CDedispSearch::m_ResultsPath.c_str(),x,y);
//    CMWAFits& dynspec_in = dynspec;

    int write_ret = dynspec.WriteFits( out_full_path );
    if( write_ret ){
        printf("ERROR : could not write dynamic spectrum FITS file %s (pixel (%d,%d))\n",out_full_path,x,y);
    }
    printf("DEBUG : analysing pixel (%d,%d) with mean=%.2f, rms=%.2f dynamic spectrum saved to file %s\n",x,y,dynspec.m_Mean,dynspec.m_RMS,out_full_path);
                
//    double threshold_in_sigma = 2;
    vector<int> good_columns;
    int pixel_count = dynspec.CutValues( threshold_in_sigma, good_columns ); // 2 sigma cut 
    int good_columns_count = good_columns.size();
    
    if ( good_columns.size() <=0 ){
       printf("WARNING : no good columns identified -> exiting now\n");
       return -1;
    }
    
// at the moment not needed, but might be used :
//    dynspec.GetStat(  dynspec.m_Mean, dynspec.m_RMS, dynspec.m_Min, dynspec.m_Max, 2, 2, good_columns[good_columns_count-1], dynspec.GetYSize(),
//                      ZERO_VALUE );

    
    // skip ZERO columns :
    sprintf(out_full_path,"%s/dynamic_spectrum_search/%04d_%04d_cut.fits",CDedispSearch::m_ResultsPath.c_str(),x,y);
    
    // TESTING saving of dynamic spectra without ZERO columns : if not working change parameter name back to dynspec (not dynspec_in)
//    CMWAFits dynspec( out_full_path, good_columns.size(), dynspec.GetYSize() );
//    dynspec_in.CopyWithoutBadColumns( dynspec, good_columns );
// TODO : just add a function to save without ZERO COLUMS 
    
    write_ret = dynspec.WriteFits( out_full_path );
    if( gDebugLocal >= 1 ){
        printf("DEBUG : written file %s with %d pixels above threshold %.2f sigma\n",out_full_path,int(good_columns.size()),threshold_in_sigma);
    }
    if( write_ret ){
        printf("ERROR : could not write dynamic spectrum FITS file %s (pixel (%d,%d))\n",out_full_path,x,y);
    }

    CMWAFits used( "used", dynspec.GetXSize(), dynspec.GetYSize() );                
//    for(int start_time=m_SkipTimeSteps;start_time<(m_MWADataCube.m_Timesteps-m_SkipTimeSteps);start_time++)
    for(int i=0;i<good_columns.size();i++)
    {
        int start_time = good_columns[i];
        
        if( start_time < m_SkipTimeSteps || start_time >= (m_MWADataCube.m_Timesteps-m_SkipTimeSteps) ){
           printf("WARNING : image at timeindex = %d skipped due to request to skip %d first and last images\n",start_time,m_SkipTimeSteps);
           continue;
        }
        
        // printf("DEBUG : processing start_time = %d\n",start_time);        
    
        for(map< string, CDispersedPixels >::iterator it = m_DispersionTemplates.begin(); it != m_DispersionTemplates.end(); ++it) {
            std::vector<CPixel>& pixels = it->second;
            double dm = atof( it->first.c_str() );
            if( gDebugLocal >= 2 ){
                printf("   DM=%s pc/cm^3 : %d pixels in dispersion measure\n",it->first.c_str(),(int)(pixels.size()));
            }
                        
            used.SetValue(0);
            double sum=0;
            int    sum_count=0;
            vector<CPixel> used_pixels;
            for(int p=0;p<pixels.size();p++){
                CPixel& pixel = pixels[p];
                           
                int time_index = start_time + pixel.time;
                int channel_index = pixel.cc; // OK ???
                int coarse_channel = m_MWADataCube.m_CoarseChannels[channel_index];
                double freq_cc_mhz = coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00;
                if( p_channels_list && p_channels_list->size()>coarse_channel ){
                   freq_cc_mhz = (*p_channels_list)[channel_index];
                }
                if( gDebugLocal >= 3 ){
                    printf("\t\tDEBUG pixel(%d) : time_index = %d, cc= %d (%d) -> freq = %.2f MHz\n",p,time_index,channel_index,coarse_channel,freq_cc_mhz);
                }
                                                     
                double max_value = -1e6;
                int max_cc=-1;
                int max_tt=-1;
                           
                for(int tt=(time_index-max_radius);tt<(time_index+max_radius);tt++){
                    for(int cc=(channel_index-max_radius);cc<(channel_index+max_radius);cc++){
                        if( tt>=0 && tt<dynspec.GetXSize() && cc>=0 && cc<dynspec.GetYSize() ){
                            if( used.getXY(tt,cc) <= 0 ){
                                double val = dynspec.getXY( tt , cc );
                                if( val > max_value && val > ZERO_VALUE ){
                                    max_value = val;
                                    max_cc    = cc;
                                    max_tt    = tt;
                                }                                      
                            }
                        }
                    }
                }
               
                if( max_cc >= 0 && max_tt >= 0 && max_value > ZERO_VALUE ){ // max_value > ZERO_VALUE means the pixel is the one above
                                                                            // the threshold used to set ZERO_VALUE - 2sigma or whatever 
                    used.setXY( max_tt, max_cc, 1 );
                    used_pixels.push_back( CPixel(max_tt,max_cc,max_value) );
                              
                    sum += max_value;
                    sum_count++;
                }
                           
                // double val = dynspec.getXY( time_index, channel_index );
                // printf("\t\tPixel[%d,%d] = %.4f\n",time_index,channel_index,val);
            }                       
                       
            if( sum_count > min_good_pixels ){
                double mean_val = (sum/sum_count);
                double rms_mean = dynspec.m_RMS / sqrt(sum_count);
                double snr = mean_val / rms_mean;
                
                if( gDebugLocal >= 3 ){ // || ( fabs(x-64)<=1 && fabs(y-64)<=1 ) ){
                   printf("DEBUG (%d,%d) : start_time=%d , dm=%.2f : sum = %.2f , sum_count = %d -> mean_value = %.2f -> snr = %.4f / %.4f = %.2f\n",x,y,start_time,dm,sum,sum_count,mean_val,mean_val,rms_mean,snr);
                }
                
                if( snr > m_SNRThreshold ){  
                   CEvent evt( start_time, x, y, eAlgoTypeTemplates, dm, snr, mean_val );
                   evt.m_Pixels = used_pixels;                          
//                 m_Events.add( start_time, x, y, eAlgoTypeTemplates, dm, snr, mean_val  );

                   if( candidates_file ){
                       mystring szUsedPixels;
                       for(int p=0;p<used_pixels.size();p++){
                           szUsedPixels << "(" << used_pixels[p].time << "," << used_pixels[p].cc << "," << used_pixels[p].flux << "),";
                       }
                   
                       candidates_file->Printf("%03d               %03d %03d   %.2f     %.2f  %.2f   %d:%s\n",evt.start_timeidx,evt.x,evt.y,evt.flux,evt.dm,evt.snr,evt.m_Pixels.size(),szUsedPixels.c_str());
                   }
                   if( gDebugLocal >= 3 ){
                      printf("%d        %d     %d       %.2f     %.2f      %.2f\n",evt.start_timeidx,evt.x,evt.y,evt.flux,evt.dm,evt.snr);
                   }
                   m_Events.push_back( evt ); 
                }
            }else{
                if( gDebugLocal >= 1 ){
                    printf("DEBUG : number of good pixels < %d\n",min_good_pixels);
                }
            }
        }       
    }
    
    return m_Events.size();
}

int CDedispSearch::FindTransientsDynSpec( double threshold_in_sigma, int border, int max_radius /* = 1 */, int min_good_pixels )
{
    double inttime=0.5;
    if( m_MWADataCube.m_Metafits ){
       inttime = m_MWADataCube.m_Metafits->inttime;
    }else{
       printf("WARNING : m_MWADataCube.m_Metafits = NULL -> cannot get real integration time !\n");
    }
    
    // quick and to match python version of the code :
    double freq_upper_end_mhz = m_MWADataCube.m_FreqUpperMHz;    
    double delta_time         = inttime;

    // loop over pixels - to check dynamic spectra of all pixels :
    
    
    // Get image at channel=0, timestep=0 so that I know (X,Y) size of the image and 
    // then I go over all the pixels and analyse their dynamic spectrum :
    // then the algorithm goes over all (X,Y) pixels in the data cube and finds transient candidates.
    // So at the end of this function (after all the loops) I have all the candidates in m_Events :    
    CMWAFits* pImage = m_MWADataCube.GetImage( m_MWADataCube.m_FirstCorrectChannel, m_MWADataCube.m_FirstCorrectTimestep );

    
    // removing all the events from m_Events list 
    m_Events.clear();
    
    if( pImage ){        
        char out_full_path[1024];
        sprintf(out_full_path,"%s/dynamic_spectrum_search/candidates_dynaspec.log",CDedispSearch::m_ResultsPath.c_str());
        string szMode = "a+";
        if( m_bForceOverwrite ){ szMode = "w"; }  
        MyOFile candidates_file( out_full_path, szMode.c_str() );
        candidates_file.Printf( LOG_FILE_HEADER );

    
        int x_size = pImage->GetXSize(), y_size = pImage->GetYSize();
        CMWAFits dynspec( "dynaspec", m_MWADataCube.m_Timesteps, m_MWADataCube.m_Channels );

        printf("DEBUG : analysing y-s in range ( %d - %d )\n",border,(y_size-border));        
        for(int y=border;y<(y_size-border);y++){
            printf("Progress : looking for transients at y = %d / %d\n",(y-border),(y_size-2*border));
            for(int x=border;x<(x_size-border);x++){
                m_MWADataCube.GetDynamicSpectrum( x, y, dynspec );                
                
                int ret = AnalyseDynSpec( dynspec, x, y, threshold_in_sigma, &candidates_file, max_radius, min_good_pixels  );
                
            }
        }
    
    
        // char out_full_path[1024];
        /*sprintf(out_full_path,"%s/dynamic_spectrum_search/candidates_dynaspec.log",CDedispSearch::m_ResultsPath.c_str());
        MyOFile candidates_file( out_full_path, "w" );
        candidates_file.Printf("# START_TIMESTEP X Y Flux[Jy]    DM    SNR\n");
    
        CMWAFits dynspec_evt( "dynaspec_candidate", m_MWADataCube.m_Timesteps, m_MWADataCube.m_Channels ); 
        for( CEventList::iterator evt=m_Events.begin(); evt != m_Events.end(); evt++ ){
            candidates_file.Printf("%d   %d   %d   %.2f     %.2f      %.2f\n",evt->start_timeidx,evt->x,evt->y,evt->flux,evt->dm,evt->snr);
            printf("%d   %d   %d   %.2f     %.2f      %.2f\n",evt->start_timeidx,evt->x,evt->y,evt->flux,evt->dm,evt->snr);
            
            // save dynamic spectrum :
            m_MWADataCube.GetDynamicSpectrum( evt->x, evt->y, dynspec_evt );

            char out_full_path[1024];
            sprintf(out_full_path,"%s/dm_%06.2f/dynamic_spectrum_search/start_timeindex%04d_%04d_%04d_%04dsigmas.fits",CDedispSearch::m_ResultsPath.c_str(),evt->dm,evt->start_timeidx,evt->x,evt->y,(int)evt->snr);
            int write_ret = dynspec.WriteFits( out_full_path );
            if( write_ret ){
                printf("ERROR : could not write dynamic spectrum FITS file %s (pixel (%d,%d))\n",out_full_path,evt->x,evt->y);
            }                          
        }
        printf("info : saved %d events to log file %s\n",m_Events.size(),out_full_path);*/
    }else{
       printf("ERROR in code : pImage = NULL in CDedispSearch::FindTransientsDynSpec\n");       
    }
    
    return m_Events.size();
}

void CDedispSearch::PrintTemplates()
{
   printf("De-dispersion templates:\n");
   for(map< string, CDispersedPixels>::iterator it = m_DispersionTemplates.begin(); it != m_DispersionTemplates.end(); ++it) {
      if( gDebugLocal >= 2 ){
          printf("   DM=%s pc/cm^3 : %d pixels in dispersion measure\n",it->first.c_str(),(int)(it->second.size()));
      }
   }   
}

int CDedispSearch::CalcTemplates( double inttime /* = -1.00 */ , double freq_upper_end_mhz /* = -1 */ , std::vector<double>* p_channels_list /*=NULL*/ ){
    return CalcTemplates( m_MWADataCube , inttime , freq_upper_end_mhz , p_channels_list );
}

int CDedispSearch::CalcTemplates( CMWAFits& dynamic_spectrum, 
                                  double inttime /* = -1 */ , std::vector<double>* p_channels_list /*=NULL*/ )
{
   printf("INFO : CDedispSearch::CalcTemplates( CMWAFits& dynamic_spectrum )\n");
   int n_channels  = dynamic_spectrum.GetYSize();
   int n_timesteps = dynamic_spectrum.GetXSize();
   
      // if inttime not specified :
   if( inttime <= 0.00 ){
       inttime = 0.5;
       if( m_MWADataCube.m_Metafits ){
           inttime = m_MWADataCube.m_Metafits->inttime;
       }
   }
//   if( freq_upper_end_mhz <= 0 ){
//       freq_upper_end_mhz = datacube.m_FreqUpperMHz;
//   }
   // quick and to match python version of the code :
   int first_coarse_channel  = m_MWADataCube.m_FirstCoarseChannel;
   double freq_lower_end_mhz = first_coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // MWA_COARSE_CHANNEL/2.00 + m_MWADataCube.m_CoarseChannels[0]*MWA_COARSE_CHANNEL;
   double freq_upper_end_mhz = freq_lower_end_mhz + (MWA_COARSE_CHANNELS-1)*MWA_COARSE_CHANNEL; // 156*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // m_MWADataCube.m_FreqUpperMHz;
   double delta_time         = inttime;
   double min_allowed_rms    = 0.00001;
   
   printf("Calculating template for DM range %.3f - %.3f and inttime = %.6f [sec] , freq_upper_end_mhz = %.3f [MHz], freq_lower_end_mhz = %.3f [MHz]\n",m_MinDM,m_MaxDM,inttime,freq_upper_end_mhz,freq_lower_end_mhz);
    
   double dm = m_MinDM;

   while( dm <= m_MaxDM ){
       char szDM[16];
       sprintf(szDM,"%03.1f",dm);
       CDispersedPixels dispersion_pixels;
       printf("Filling template for DM = %s ( n_channels = %d  )\n",szDM,n_channels);
  
       for(int start_timeindex=0;start_timeindex<n_timesteps;start_timeindex++){
           double unixtime = m_MWADataCube.m_StartUnixTime + start_timeindex * inttime;
        
//        CMWAFits* pStartImage = m_MWADataCube.GetImage( 0, start_timeindex );
//        if( !pStartImage ){
//           continue;
//        }
           double start_unixtime = 0.00; // m_MWADataCube.m_StartUnixTime;
//        if( pStartImage ){
//           start_unixtime = pStartImage->GetUnixTime();
//        }
         
           int used_channels = 0;        
           double arrival_time = 0.00 , arrival_time_step = 0.50;
           double freq_cc_mhz = freq_upper_end_mhz;
           int timeindex = -1;
           double min_dist = 0.00;
        
           // this loop goes over arrival times 
           while( freq_cc_mhz >= freq_lower_end_mhz && arrival_time<=(inttime*n_timesteps) && timeindex < n_timesteps ){
               double freq_mhz = freq_cc_mhz;
               int coarse_channel = int( freq_mhz / MWA_COARSE_CHANNEL ); // floor -> absolute number of coarse channel 
               int cc = ( coarse_channel - first_coarse_channel);
               
               // find index in the list of coarse channels in the metadata :
//            int cc = find_value( m_MWADataCube.m_CoarseChannels, coarse_channel ); // coarse_channel - m_MWADataCube.m_CoarseChannels[0]; // WARNING : assumes it's not a picket fence observation !
            
               if( cc < 0 ){
                  printf("ERROR : freq_mhz = %.2f MHz -> coarse channel = %d not found in the observation metadata - this means bug in the code (as it should not happen) -> cannot continue / exiting function\n",freq_mhz,coarse_channel);
                  return -1;
               }
            
               if( cc < 24 ){        
//                   double unixtime_fits = GetUnixTime();                
                   if( gDebugLocal >= 2 ){
                       printf("\tCalculating arrival time from upper end of frequency at %.4f MHz to upper end of channel = %d ( %d / %.4f MHz )\n",freq_upper_end_mhz,cc,coarse_channel,freq_cc_mhz);
                   }
            
                   double dispersion_delay_ms = arrival_time*1000; //  dispersion_milisec( freq_upper_end_mhz, freq_cc_mhz, dm ); // # freq_mhz+23.00*1.28+0.5*1.28 - is the higher end of freq. band 
                   double dispersion_delay_sec = dispersion_delay_ms / 1000.00;
                   double dispersion_delay_timesteps_float = dispersion_delay_sec / delta_time;

// unitl I started to implement multi-obsID version it was int() , but decided to start testing round() :
//                int dispersion_delay_timesteps_int = int( dispersion_delay_timesteps_float );
// Nope I think it was int() on purpose, because we want to include image when the burst started not the next one as in that case it would only 
// containt signal if considerably scattered. 
// SO : it should indeed be int()
// But how to implement int() in m_MWADataCube.GetClosestImage ???
                   int dispersion_delay_timesteps_int = int( dispersion_delay_timesteps_float );
            
                   // time index :
                   int timeindex_estimate = start_timeindex + dispersion_delay_timesteps_int;
                   // Proper calculation using time :
                   double image_unixtime = start_unixtime + dispersion_delay_sec;
                   timeindex = m_MWADataCube.GetOptimalImage( image_unixtime , min_dist ); // was GetClosestImage
                   if( timeindex <= 0 ){
                      timeindex = timeindex_estimate;
                      min_dist  = 0.00;
                   }
// WORKS OK :                int timeindex = timeindex_estimate;
                   if( gDebugLocal >= 2 ){
                      printf("DEBUG : image_unixtime = %.2f ( = %.2f + %.2f ) -> closest timeindex = %d ( estimated guess = %d = %d + %d ) , min_dist = %.2f [sec]\n",image_unixtime,m_MWADataCube.m_StartUnixTime,dispersion_delay_sec,timeindex,timeindex_estimate,start_timeindex , dispersion_delay_timesteps_int,min_dist);
                   }
            
                   if( timeindex < n_timesteps && timeindex >= 0 && min_dist <= inttime ){
                      CPixel pixel;
                      pixel.time = dispersion_delay_timesteps_int;
                      pixel.cc   = cc;
                      dispersion_pixels.push_back( pixel );
                   }else{                    
                       if( gDebugLocal>=2){
                          printf("\tSKIPPED : cc = %.2f MHz : delay = %.2f [sec] -> delta_timesteps=%.2f=%d >= %d timesteps , timeindex = %d , min_dist = %.4f [sec]",freq_cc_mhz,dispersion_delay_sec,dispersion_delay_timesteps_float,dispersion_delay_timesteps_int,m_MWADataCube.m_Timesteps,timeindex,min_dist);                       
                       }
                       if( timeindex < 0 ){
                          printf("\t\tWARNING : timeindex = %d !!!\n",timeindex);
                       }
                   }
               }else{
                  if( gDebugLocal>=2){printf("WARNING : missing image at coarse channel = %d\n",cc);}
               }
            
               arrival_time += arrival_time_step;
               freq_cc_mhz = dispersion_freq_mhz( freq_upper_end_mhz , arrival_time*1000.00, dm );
           }        
       }    
       
       m_DispersionTemplates[szDM] = dispersion_pixels;
       printf("\tSet pixels for DM = %s\n",szDM);
       dm = dm + m_StepDM;
   }
   
   return 1;
}


int CDedispSearch::CalcTemplatesExact( CMWAFits& dynamic_spectrum )
{
   printf("INFO : CDedispSearch::CalcTemplatesExact( CMWAFits& dynamic_spectrum )\n");
   if( !m_MWADataCube.m_Metafits ){
       printf("ERROR : to correctly generate time series images metafits file is required !\n");
       exit(-1);
   }
   
   // remove old :
   m_DispersionTemplates.clear();
   printf("Cleared DM templated , current size = %d\n",int(m_DispersionTemplates.size()));

 
   int n_channels  = dynamic_spectrum.GetYSize();
   int n_timesteps = dynamic_spectrum.GetXSize();
   double inttime = m_MWADataCube.m_Metafits->inttime;
   int n_coarse_channels = m_MWADataCube.m_CoarseChannels.size();
   double first_coarse_channel  = m_MWADataCube.m_CoarseChannels[0];
   double last_coarse_channel   = m_MWADataCube.m_CoarseChannels[n_coarse_channels-1];
//   double freq_lower_end_mhz = first_coarse_channel*MWA_COARSE_CHANNEL - MWA_COARSE_CHANNEL/2.00; // MWA_COARSE_CHANNEL/2.00 + m_MWADataCube.m_CoarseChannels[0]*MWA_COARSE_CHANNEL;
//   double freq_upper_end_mhz = last_coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // m_MWADataCube.m_FreqUpperMHz;
   double freq_lower_end_mhz = m_MWADataCube.m_fCoarseChannelFreqMHz[0] - MWA_COARSE_CHANNEL/2.00; 
   double freq_upper_end_mhz = m_MWADataCube.m_fCoarseChannelFreqMHz[n_coarse_channels-1] + MWA_COARSE_CHANNEL/2.00;
   
   printf("DEBUG : n_channels = %d , n_timesteps = %d , inttime = %.8f [sec] , n_coarse_channels = %d , first_coarse_channel = %.4f , freq_lower_end_mhz = %.8f [MHz] , last_coarse_channel = %.4f , freq_upper_end_mhz = %.8f [MHz]\n",
                 n_channels,n_timesteps,inttime,n_coarse_channels,first_coarse_channel,freq_lower_end_mhz,last_coarse_channel,freq_upper_end_mhz);
   
   if( n_channels != n_coarse_channels ){
      printf("ERROR : dynamic spectrum has %d channels, but meta data has %d coarse channels -> cannot continue\n",n_channels,n_coarse_channels);
      return -1;
   }
   
//   double delta_time = 0.001;   

//   if( freq_upper_end_mhz <= 0 ){
//       freq_upper_end_mhz = datacube.m_FreqUpperMHz;
//   }
   // quick and to match python version of the code :
//   int first_coarse_channel  = 133;
//   double freq_lower_end_mhz = first_coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // MWA_COARSE_CHANNEL/2.00 + m_MWADataCube.m_CoarseChannels[0]*MWA_COARSE_CHANNEL;
//   double freq_upper_end_mhz = 156*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; // m_MWADataCube.m_FreqUpperMHz;
//   double delta_time         = inttime;
   double min_allowed_rms    = 0.00001;
   double max_sweep_time = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, m_MaxDM ) / 1000.00; // in seconds
   
   printf("Calculating template for DM range %.3f - %.3f and inttime = %.6f [sec] , freq_upper_end_mhz = %.6f [MHz], max_sweep_time = %.2f [ms] , coarse channels %d x %.2f MHz\n",m_MinDM,m_MaxDM,inttime,freq_upper_end_mhz,(max_sweep_time*1000.00),MWA_COARSE_CHANNELS,MWA_COARSE_CHANNEL);
    
   // 2023-07-17 - in fact start DM should be centre of first DM bin , so for example 0 - 150 -> 75 !
   // double dm = m_MinDM;
   // also requires change in ::GetDedispersedTimeSeries : int out_dm_index = round( (dm-m_MinDM-m_StepDM/2.00) / m_StepDM);
   double dm = m_MinDM + m_StepDM/2.00; 
   
   while( dm <= m_MaxDM ){
       bool bDebug=false;
       char szDM[16];
       sprintf(szDM,"%06.1f",dm);
       CDispersedPixels dispersion_pixels, dispersion_pixels_all;
       printf("Filling template for DM = %s ( n_channels = %d  , %.6f - %.6f MHz )\n",szDM,n_channels,freq_upper_end_mhz,freq_lower_end_mhz);
       
       double sweep_time_dm = dispersion_milisec( freq_upper_end_mhz, freq_lower_end_mhz, dm )/1000.00; // in seconds
       double delta_time = 0.00001; // second vs. 0.001 if too large for low-DM objects steps are too large and not all channles are added !!!
       
//       double arrival_time = 0.00 , arrival_time_step = 0.50;
//       double min_dist = 0.00;

       if( fabs(dm-2.9) <= 0.00001 ){
          printf("DEBUG : of dispersion paths enabled\n");          
          bDebug = true;
       }
       
       
       // go through time in small steps to find all the pixels belonging to DM sweep :
       // WARNING : it does not work for DM = 0 I am just adding 593 (all pixels along time axis !!!)
       double unixtime = 0.00;
       int timeindex = -1;
       int used_channels = 0;         
       double freq_mhz = freq_upper_end_mhz;
       double prev_freq_mhz = freq_upper_end_mhz;
       double freq_center_upper_channel_mhz = m_MWADataCube.m_CoarseChannels[n_coarse_channels-1]*MWA_COARSE_CHANNEL;
       double lowest_coarse_channel = m_MWADataCube.m_CoarseChannels[0];
       int last_added_timeindex = -1;
       int last_added_cc = 1e9;
       
       if( dm < 0.000000000001 ){ // == 0 
          // dm = 0 -> just add all channels :
          for(int cc=(n_channels-1);cc>=0;cc--){
             int diff_cc = (n_channels-1) - cc;
             CPixel pixel;
             pixel.time = unixtime;
             pixel.cc = cc;
             pixel.exact_uxtime = unixtime;
             pixel.exact_freq_mhz = m_MWADataCube.m_fCoarseChannelFreqMHz[n_coarse_channels-1] - diff_cc*MWA_COARSE_CHANNEL;
             dispersion_pixels.push_back( pixel );
             dispersion_pixels_all.push_back( pixel );
             printf("DEBUG : DM = %.4f -> added pixel (%d,%d), freq = %.4f MHz , current count = %d\n",dm,pixel.time,pixel.cc,pixel.exact_freq_mhz,int(dispersion_pixels.size()));
          }
          
       }else{                      
          while( unixtime <= n_timesteps*inttime && freq_mhz >= freq_lower_end_mhz && timeindex < n_timesteps && (last_added_cc - lowest_coarse_channel) >=0 ){
             double delta_time_total_ms = ( unixtime )*1000.00;
             double freq_arrival_end_mhz = m_MWADataCube.m_fCoarseChannelFreqMHz[n_coarse_channels-1]; // use centre of the highest channel not the edge ( was freq_upper_end_mhz );
             freq_mhz = dispersion_freq_mhz( freq_arrival_end_mhz , delta_time_total_ms, dm ); // 2023-05-25 : was freq_center_upper_channel_mhz 
             double delta_freq_mhz = freq_mhz - prev_freq_mhz;

             if( unixtime >= 0 ){
                // We only care about times >=0 to pick up pixels            
                // 2023-05-25 : if we assume very narrow pulse (inifitly narrow) than it should rather be int() than round()
                timeindex = int( unixtime / inttime ); //2021-04-09 : was 0.5 -> inttime      // rounded down because time starts from 0.00 and that's when we enter first bin time, we need 0.5 seconds to get to next bin 
                                                      // 2023-02-22 : int -> round 
                                                      // 2023-05-25 : temporary back from round -> int as the bin is really where the signal is ...
                int cc        = round( freq_mhz /  MWA_COARSE_CHANNEL ); // 2021-04-09 : was 1.28 -> MWA_COARSE_CHANNEL to nearest integer , 
                                                                      // 2023-02-22 : int -> round
                                                                      // 2023-05-25 : temporary back from round -> int as the bin is really where the signal is ...
                                                                      //            !!! round because channels ch=121 -> center frequency is 121*154.88 
             
                if( bDebug ){
                   printf("DEBUG : dm = %.3f : timeindex = %d, cc = %d ( unixtime = %.8f , freq_mhz = %.4f MHz )\n",dm,timeindex,cc,unixtime,freq_mhz);
                }
             
                CPixel pixel;    
                if( cc != last_added_cc || timeindex != last_added_timeindex ){
                   pixel.time = timeindex;
                   pixel.cc   = cc - lowest_coarse_channel; // adding relative to first coarse channel 
                
                   if( pixel.cc < 0 || pixel.cc >= n_channels ){
                      printf("WARNING : calculated channel outside 0 - %d range\n",n_channels);
                   }

/*                if( pixel.cc < 0 ){
                   printf("ERROR in code : coarse channel cannot be negative !!! cc = %d\n",pixel.cc);
                }*/
                                
                   if( pixel.cc >= 0 ){
                      // we only add when above lower frequency edge :
                      dispersion_pixels.add_unique( pixel );
                      printf("DEBUG : DM = %.4f -> added pixel (%d,%d), freq = %.4f MHz, current count = %d\n",dm,pixel.time,pixel.cc,freq_mhz,int(dispersion_pixels.size()));
                   }
                
                   last_added_timeindex = timeindex;
                   last_added_cc = cc;
                }

             
                // fill list of all pixels :
                pixel.time = timeindex;
                pixel.cc   = cc - lowest_coarse_channel; // adding relative to first coarse channel           
                pixel.exact_uxtime = unixtime;
                pixel.exact_freq_mhz = freq_mhz;             
                dispersion_pixels_all.push_back( pixel );
             }
              
             prev_freq_mhz = freq_mhz;           
             unixtime += delta_time;
          }                 
       }
       
       // TODO : to debug print exact_uxtime and exact_freq_mhz from dispersion_pixels_all to then plot the path and see 
       // how to calculate weights.
       // 2023-05-22 :
       // Create 2D array with the same number of channels and timestamps as the dynamic spectrum and count "hits" of exact_uxtime and exact_freq_mhz per pixel 
       // then calculate weights for time bins within a single freq. channel -> will be correct !!!
       // apply these weigths to the pixels in the array dispersion_pixels
       CBgFits counter( last_added_timeindex+10 , n_channels );
       CBgFits weighted( last_added_timeindex+10 , n_channels );
       CBgFits generated( last_added_timeindex+10 , n_channels );
       counter.SetValue(0.00);
       
       for(int pix=0;pix<dispersion_pixels_all.size();pix++){
          // calculate time bin index and freq. channel 
          CPixel& pixel = dispersion_pixels_all[pix];
          
          int time_bin_index = int( pixel.exact_uxtime / inttime ); // int or round - currently consistent with the earlier code above to calculate pixel.cc and pixel.time 
          double diff_freq = freq_upper_end_mhz - pixel.exact_freq_mhz;
          
          int ch = round( pixel.exact_freq_mhz / MWA_COARSE_CHANNEL ); // calculate which freq. channel bin it falls into. int or round - currently consistent with the earlier code above to calculate pixel.cc and pixel.time
          int ch_diff = ch - lowest_coarse_channel;
          
          if( CDedispSearch::m_bDebugDedispersion ){
             printf("DEBUG : (%.3f,%.3f) / %.3f / %.3f = %d , %d\n",pixel.exact_uxtime,pixel.exact_freq_mhz,inttime,MWA_COARSE_CHANNEL,time_bin_index,ch);
          }
          counter.addXY( time_bin_index, ch_diff, 1 );
       }  
       counter.WriteFits("counter.fits");

       // calculate weights:       
       for(int ch=0;ch<n_channels;ch++){
          double total_counts_ch = counter.Sum(ch);
          
          counter.Divide( ch, total_counts_ch );
       }
       counter.WriteFits("weights.fits");
       
       // calculate verification image :
       generated.SetValue(0.00);
       for(int pix=0;pix<dispersion_pixels.size();pix++){
          CPixel& pixel = dispersion_pixels[pix];
          
          if( pixel.cc >=0 && pixel.cc < generated.GetYSize() && pixel.time>=0 && pixel.time < generated.GetXSize() ){
             generated.addXY( pixel.time, pixel.cc, 1.00 );
          }
       }
       generated.WriteFits("template.fits");
       
       
       
       // weight pixels by weights :
       // dispersion_pixels.
       weighted.SetValue(0.00);
       for(int pix=0;pix<dispersion_pixels.size();pix++){
          CPixel& pixel = dispersion_pixels[pix];
          double pixel_weight = counter.getXY( pixel.time, pixel.cc );
          pixel.weight = pixel_weight;

          if( pixel.cc >=0 && pixel.cc < counter.GetYSize() && pixel.time>=0 && pixel.time < counter.GetXSize() ){
             weighted.addXY( pixel.time, pixel.cc, 1.00*pixel_weight );
          }
       }
       counter.WriteFits("template_weighted.fits");

       
       m_DispersionTemplates[szDM] = dispersion_pixels;
       
       // saving 
       char reg_filename[512];
       sprintf( reg_filename , "dispersed_path_dm%.3f.reg", dm );
       dispersion_pixels.save_reg_file( reg_filename );      
       
       sprintf( reg_filename , "dispersed_path_dm%.3f_all.txt", dm );
       dispersion_pixels_all.save_text_file( reg_filename );
       
       printf("\tSet pixels for DM = %s ( %d dispersion pixels)\n",szDM,int(dispersion_pixels.size()));
       dm = dm + m_StepDM;
   }
   
   return 1;
}


int CDedispSearch::CalcTemplates( CMWADataCube& datacube, double inttime /* = -1 */ , double freq_upper_end_mhz /* = -1 */ , std::vector<double>* p_channels_list /*=NULL*/ )
{
   printf("INFO : CDedispSearch::CalcTemplates( CMWADataCube& datacube )\n");
   int n_channels = datacube.m_Channels;
   int n_timesteps = datacube.m_Timesteps;
   
   // if inttime not specified :
   if( inttime <= 0.00 ){
       inttime = 0.5;
       if( datacube.m_Metafits ){
           inttime = datacube.m_Metafits->inttime;
       }
   }
   if( freq_upper_end_mhz <= 0 ){
       freq_upper_end_mhz = datacube.m_FreqUpperMHz;
   }
   printf("Calculating template for DM range %.3f - %.3f and inttime = %.6f [sec] , freq_upper_end_mhz = %.3f [MHz]\n",m_MinDM,m_MaxDM,inttime,freq_upper_end_mhz);

   int start_timeindex=0;
   double unixtime = datacube.m_StartUnixTime + start_timeindex * inttime;
  
   double dm = m_MinDM;
   
   while( dm <= m_MaxDM ){
       char szDM[16];
       sprintf(szDM,"%03.1f",dm);
       CDispersedPixels dispersion_pixels;
       printf("Filling template for DM = %s ( n_channels = %d  )\n",szDM,n_channels);
   
       for(int cc=0;cc<n_channels;cc++){
           if( cc < m_MWADataCube.m_CoarseChannels.size() ){
               int coarse_channel = datacube.m_CoarseChannels[cc];
               double freq_cc_mhz = coarse_channel*MWA_COARSE_CHANNEL + MWA_COARSE_CHANNEL/2.00; //  should be added to calculate arrival time at the
               if( p_channels_list && p_channels_list->size()>=n_channels ){
                  freq_cc_mhz = (*p_channels_list)[cc];
               }
                
               double dispersion_delay_ms = dispersion_milisec( freq_upper_end_mhz, freq_cc_mhz, dm ); // # freq_mhz+23.00*1.28+0.5*1.28 - is the higher end of freq. band 
               double dispersion_delay_sec = dispersion_delay_ms / 1000.00;
               double dispersion_delay_timesteps_float = dispersion_delay_sec / inttime;
               int dispersion_delay_timesteps_int = int( dispersion_delay_timesteps_float );
               
               printf("   cc=%d -> freq = %.2f MHz -> delay = %.2f ms = %.2f steps\n",cc,freq_cc_mhz,dispersion_delay_ms,dispersion_delay_timesteps_float);
               
               CPixel pixel;
               pixel.time = dispersion_delay_timesteps_int;
               pixel.cc   = cc;
               dispersion_pixels.push_back( pixel );
            }
        }
        
        m_DispersionTemplates[szDM] = dispersion_pixels;
        printf("\tSet pixels for DM = %s\n",szDM);
        dm = dm + m_StepDM;
   }
   
   return 1;
}

int CDedispSearch::FindTransients( double dm, double threshold_in_sigma, int border )
{
    int n_candidates = 0;
    CBgFits dynspec( m_MWADataCube.m_Timesteps, m_MWADataCube.m_Channels );
    
    printf("CDedispSearch::FindTransients : dm=%.2f , threshold_in_sigma=%.2f , border=%d\n",dm, threshold_in_sigma, border);    

    char out_full_path[1024];
    sprintf(out_full_path,"%s/dm_%06.2f/candidates_dm%05.2f.log",CDedispSearch::m_ResultsPath.c_str(),dm,dm);
    string szMode = "a+";
    if( m_bForceOverwrite ){ szMode = "w"; }  
    bool bPrintHeader = true;
    if( MyFile::DoesFileExist( out_full_path ) ){
       bPrintHeader = false;
    }
    MyOFile candidates_file( out_full_path, szMode.c_str() );
    if( bPrintHeader ){
        candidates_file.Printf("# START_TIMESTEP X Y Flux[Jy]    DM   SNR\n");
    }
//    WARNING : when using a single CBgFits object to save all the files 
//              it crashes when m_fptr != NULL (so either uncomment line //      m_fptr = NULL; // NEW 2016-09-28 - )
//              or somehow set m_fptr to NULL 
    
    for(int i=0;i<m_DedispersedImages.size();i++){
        CMWAFits* pDedispImage = m_DedispersedImages[i];
        
 
        if( pDedispImage ){
            if( pDedispImage->m_UsedChannels > CDedispSearch::m_MinUsedChannels && (pDedispImage->m_OriginalImages.size()<=0 || pDedispImage->m_OriginalImages.size() > m_MinUsedImagesInDedisp) ){
                double mean = pDedispImage->m_Mean;
                double rms   = pDedispImage->m_RMS;
                double threshold_jy = threshold_in_sigma * rms + mean;
            
                int x_size = pDedispImage->GetXSize();
                int y_size = pDedispImage->GetYSize();
            
                int count_above = 0;
                int accepted    = 0;
                int rejected_due_to_ref_sources = 0;
                for(int y=border;y<(y_size-border);y++){
                   for(int x=border;x<(x_size-border);x++){               
                       double val = pDedispImage->getXY(x,y);
                       bool bShowCandidate=false;
                       int event_index = m_MWADataCube.m_StartTimeIndex + i;
                   
                       if( val > threshold_jy ){                   
                           double val_in_sigmas = (val/rms);
                           double debug_value = pDedispImage->getXY(0,0);
                           count_above++;
                           if( gDebugLocal >= 1 || fabs(dm-50.00)<0.1 ){ // save all or just 1 selected
                               printf("Candidate at (%d,%d), i=%d , dm=%.2f, debug_value=%.4f, flux = %.2f Jy ( %.2f sigma ) -> ",x,y,event_index,dm,debug_value,val,val_in_sigmas);
                               bShowCandidate = true;
                           }
                       
                           CSource* pRefSource = NULL;
                           if( m_MWADataCube.m_pMeanImage ){
                               pRefSource = (m_MWADataCube.m_pMeanImage->m_Sources).find( x, y, 5 );
                           }
                           if( !pRefSource && m_MWADataCube.m_pMedianImage ){
                              pRefSource = (m_MWADataCube.m_pMedianImage->m_Sources).find( x, y, 5 );
                           }
                       
                           if( !pRefSource || m_bRejectRefSources<=0 ){
                              n_candidates++;
                              if( gDebugLocal >= 1 || bShowCandidate ){
                                  printf("Accepted\n");
                              }
                              string szComment;
                              if( pRefSource ){
                                 szComment = "REF-SOURCE";
                              }
                              
                              candidates_file.Printf("%d %d %d %.2f %.2f %.2f %s\n",event_index,x,y,val,dm,val_in_sigmas,szComment.c_str());
                          
                              // saving dynamic spectrum                          
                              char out_full_path[1024];
                              sprintf(out_full_path,"%s/dm_%06.2f/dynamic_spectrum/start_timeindex%04d_%04d_%04d_%04dsigmas.fits",CDedispSearch::m_ResultsPath.c_str(),dm,event_index,x,y,(int)val_in_sigmas);
                              m_MWADataCube.GetDynamicSpectrum( x, y, dynspec );
                              int write_ret = dynspec.WriteFits( out_full_path );
                              if( write_ret ){
                                 printf("ERROR : could not write dynamic spectrum FITS file %s (pixel (%d,%d))\n",out_full_path,x,y);
                              }
                              accepted++;
                           }else{
                              if( gDebugLocal >= 1 || bShowCandidate ){
                                  printf("Skipped due to reference source at (%.1f,%.1f) , flux = %.2f Jy\n",pRefSource->x,pRefSource->y,pRefSource->flux);
                              }
                              rejected_due_to_ref_sources++;
                           }
                       }
                   }
                }            
                printf("CDedispSearch::FindTransients : number of pixels above threshold = %.2f Jy is %d ( rej-due-refsources = %d , accepted = %d )\n",threshold_jy,count_above,rejected_due_to_ref_sources,accepted);
            }else{
                printf("WARNING : image %d (%s) skipped due to %d only channels used in de-dispersion (< minimum = %d) , or number of averaged images = %d < minimum = %d\n",i,pDedispImage->GetFileName(),pDedispImage->m_UsedChannels,CDedispSearch::m_MinUsedChannels,int(pDedispImage->m_OriginalImages.size()),m_MinUsedImagesInDedisp);
            }
        }else{
            printf("ERROR : pDedispImage = NULL - error in code ???\n");
        }        
    }
    
    return n_candidates;
}

int CDedispSearch::Run( double threshold_in_sigma, int border, int min_good_pixels )
{
   double dm = m_MinDM;

   // calculate reference images and list of sources :
   CalcReferenceImages();
   
   if( CDedispSearch::m_bRunBruteForceAlgo ){ 
      time_t start_ux = get_dttm();
      printf("Starting brute force algo at : %d\n",(int)start_ux);
      while( dm <= m_MaxDM ){       
          printf("------------------------------------------------------------ Testing dm = %.2f [pc/cm^3] ------------------------------------------------------------\n",dm);
          Dedisperse( dm );
       
          FindTransients( dm, threshold_in_sigma, border );
   
          dm = dm + m_StepDM;
      }
      
      time_t end_ux = get_dttm();
      printf("Brute force algo finished at %d (took %d sec)\n",(int)end_ux,(int)(end_ux-start_ux));
   }
   
   if( CDedispSearch::m_bRunDynSpecAlgo ){
      time_t start_ux = get_dttm();
      printf("Starting calculation of templates at : %d\n",(int)start_ux);
      // calculate templates :
      CalcTemplates();
      PrintTemplates();
      time_t end_ux = get_dttm();
      printf("Calculating templates for dynaspec algo took %d sec\n",(int)(end_ux-start_ux));

      // different way of looking for dispersed pulses :
      start_ux = get_dttm();
      printf("DEBUG : starting FindTransientsDynSpec at %d, using cut threshold = %.2f\n",(int)start_ux,CDedispSearch::m_ThresholdToCutInSigmas);fflush(stdout);
      FindTransientsDynSpec( CDedispSearch::m_ThresholdToCutInSigmas, border, 1, min_good_pixels );
      end_ux = get_dttm();
      printf("Dynaspec algo finished at %d (took %d sec)\n",(int)end_ux,(int)(end_ux-start_ux));
   }
   
   return 1;
}

