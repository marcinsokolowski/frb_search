#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include <myfile.h>
#include <mystrtable.h>
#include <myparser.h>
#include <math.h>
// #include <mystring.h>
#include "mwa_fits.h"
#include "pixel.h"
#include "dedisp_search.h"

// filterbank files :
// #include <SigprocFile.h>
#ifdef _USE_FILFILE_
#include <filfile.h>
#endif

#include <bits/stdc++.h> 
#include <vector>

int gDebugLocal=0;
int CMWAFits::m_Border=10;
bool CMWAFits::m_bUseLocalThreshold=true;

double MWA_COARSE_CHANNEL  = 1.28; // # MHz 
int    MWA_COARSE_CHANNELS = 24;  // # Number of MWA coarse channels 

double CMWAFits::GetStat( double& mean, double& rms, double& minval, double& maxval, CPixel& max_pixel, CPixel& min_pixel,
                          CBgArray& max_along_x, vector<int>& max_along_x_x,
                          CBgArray& max_along_y, vector<int>& max_along_y_y,
                          int x_start, int y_start, int x_end, int y_end,
                          double min_allowed_value, double max_allowed_value  )
{
   double sum = 0.00;
   double sum2 = 0.00;
   int    cnt  = 0;
   
   minval = 1e6;
   maxval = -1e6;
   
   if( y_end < 0 ){
      y_end = m_SizeY;
   }
   if( x_end < 0 ){
      x_end = m_SizeX;
   }
   
   max_along_x.clear();
   max_along_x.assign( m_SizeY, 0.00 );
   
   max_along_x_x.clear();
   max_along_x_x.assign( m_SizeY, 0 );

   int nan_count = 0, total_count = 0;   
   for(int y=y_start;y<y_end;y++){
    
       double max_in_line=-1e20;
       int    max_in_line_x=-1000;
       for(int x=x_start;x<x_end;x++){
           double val = valXY(x,y);
           
           total_count++;                      
           if ( isnan(val) ){
              nan_count++;
              continue;
           }

           
           if( val >= min_allowed_value && val <= max_allowed_value ){
              sum  += val;
              sum2 += val*val;
              cnt  += 1;
           
              if( val > maxval ){
                  maxval = val;
                  max_pixel.time = x;
                  max_pixel.cc = y;                  
                  max_pixel.flux = val;
              }
              if( val < minval ){
                  minval = val;
                  min_pixel.time = x;
                  min_pixel.cc = y;     
                  min_pixel.flux = val;       
              }
           }
           
           if( val > max_in_line ){
              if( x < (m_SizeX-200) ){
                 // ignore last 200 columns 
                 max_in_line = val;
                 max_in_line_x = x;              
              }
           }
       }
       
       max_along_x[y] = max_in_line;
       max_along_x_x[y] = max_in_line_x;
   }
  
   mean = 0.00;
   rms  = 0.00;
 
   if( cnt >= 1 ){
      mean = sum / cnt;
      rms  = sqrt( sum2/cnt - mean*mean );
   }
   
   if( nan_count > 0 ){
      printf("WARNING : %d / %d are NaN values found and ignored\n",nan_count,total_count);
   }

   return 0.00;
}

double CMWAFits::GetStat( int border /* = 20 */ )
{
    double ret = CBgFits::GetStatBorder( m_Mean, m_RMS, m_Min, m_Max, border );
    printf("INFO : GetStat mean = %.4f , rms = %.4f , min = %.4f , max = %.4f , border = %d\n",m_Mean, m_RMS, m_Min, m_Max, border );
    return ret;
}

// RMS per count is calculated as RMS value corresponding to number of pixels in DM-sweep
int CMWAFits::FindSources( double threshold_in_sigma, int border, const char* szOutRegFile, bool bFluxInSigmas,  int minX, int maxX, 
                           CBgFits* pCountMap, 
                           CBgArray* p_rms_per_count,
                           int max_count_param /*=-1*/
                         )
{
   printf("PROGRESS CMWAFits::FindSources size %d x %d, PARAMS : max_count_param = %d\n",int(GetXSize()),int(GetYSize()),max_count_param);
  
   GetStat( border );

   // calculate MEAN and RMS for a specific DM (vertical axis) :   
   CBgArray mean_for_dm,rms_for_dm;
   MeanLines( mean_for_dm, rms_for_dm );

   string szOutMeanForDmFile;
   change_ext( szOutRegFile, "mean_spectrum", szOutMeanForDmFile );
   mean_for_dm.SaveToFile( szOutMeanForDmFile.c_str(), NULL, &rms_for_dm );
   
   
   // RMS per number of pixels (only if pCountMap != NULL ):
   if( p_rms_per_count ){
      p_rms_per_count->clear();
   }
   
   
   if( pCountMap ){
      if( m_SizeX != pCountMap->GetXSize() || m_SizeY != pCountMap->GetYSize() ){
         printf("ERROR : size of the de-dispersed series and count pixel map have to be the same they are not (%d,%d) != (%d,%d)\n",int(m_SizeX),int(m_SizeY),int(pCountMap->GetXSize()),int(pCountMap->GetYSize()));
         exit(-1);
      }

      int max_count = m_SizeY + m_SizeX;// TODO : change to max disp-delay otherwise takes too long !!!
      if( max_count_param > 0 ){
         max_count = max_count_param;
      }
      printf("Calculating RMS per number of pixels in the map, max_count = %d ( = %d + %d )\n",max_count,int(m_SizeX),int(m_SizeY));fflush(stdout);
      int last_zeros = 0;
      int init_rms_count = max_count+5;
      if( p_rms_per_count ){
         p_rms_per_count->assign( init_rms_count, 0 );
      }
      
      // for every number of points in a sweep (count) RMS is calculated
      double max_rms=-1e20;
      for(int cnt=1;cnt<init_rms_count;cnt++){
         printf("\tcnt = %d ",cnt);fflush(stdout);
         double sum=0.00,sum2=0.00;
         int sum_count=0;
         // calculate RMS of the pixels in the de-dispersed time series based on cnt pixels in dispersion sweep :
         
         // quite expensive going through the entire DTS again to find pixels with the same number of pixels in the DM-sweep (i.e. count)
         for(int y=0;y<m_SizeY;y++){
            for(int x=0;x<m_SizeX;x++){
               int map_count = pCountMap->getXY( x, y );
               if( map_count == cnt ){
                  double fits_value = getXY( x , y );
                  sum += fits_value;
                  sum2 += fits_value*fits_value;
                  sum_count++;
               }
            }
         }
         
         double mean = 0.00, rms = 0.00;
         if( sum_count > 0 ){
            mean = sum / sum_count;
            rms = sqrt( (sum2/sum_count) - (mean*mean) );
            if( p_rms_per_count ){
               (*p_rms_per_count)[cnt] = rms;
            }
            
            if( rms > max_rms ){
               max_rms = rms;
            }
         }else{
            last_zeros++;
         }
         printf("-> rms(%d) = %.4f ( mean = %.4f , sum2 = %.4f , sum_count = %d )\n",cnt,rms,mean,sum2,sum_count);fflush(stdout);
         
         if( sum_count == 0 && last_zeros > 100 ){
            printf("DEBUG : %d last counts have count=0 -> no need to continue\n",last_zeros);
            break;
         }
         
/* not ok for many channels and large DMs - because exists too early !        
         if( sum_count == 0 && last_zeros > 50 ){ // if last 50 are zeros then stop now 
             printf("Auto-detected maximum number of pixels in path at %d -> exiting loop earlier (nominally up to %d)\n",cnt,max_count);                          
             break;
         }*/
      }
      
      // TODO :
      // implement least squars fit of PARAM0 / SQRT(N_pixels) ( starting from pixel = 3 or 4 )
      // until fit is implemented overwrite the first 2 values with the maximum value 
      if( p_rms_per_count ){
         // setting maximum RMS to 0 and 1 count (DM-sweep with 0 or 1 pixels only):
         (*p_rms_per_count)[0] = max_rms;
         (*p_rms_per_count)[1] = max_rms;      
         string szOutRmsPerCountFile;
         change_ext( szOutRegFile, "rms_per_count", szOutRmsPerCountFile );
         p_rms_per_count->SaveToFile( szOutRmsPerCountFile.c_str() );
         printf("Saved RMS per count to file %s\n",szOutRmsPerCountFile.c_str() );
      }
   }
   
   m_Sources.clear();

   if( maxX <= 0 ){
      maxX = m_SizeX;
   }
   
   int debug_maxX = -1;
   double threshold_jy_global = threshold_in_sigma*m_RMS + m_Mean; // global threshold 
   double threshold_jy = threshold_jy_global;
   printf("DEBUG GLOBAL THRESHOLD = %.4f * %.4f + %.4f = %.4f ( m_bUseLocalThreshold = %d ) -> searching for pixels in the range %d - %d\n",threshold_in_sigma,m_RMS,m_Mean,threshold_jy,m_bUseLocalThreshold,minX,maxX);
   for(int y=0;y<m_SizeY;y++){
      threshold_jy = threshold_jy_global;
      double rms_local = m_RMS;
      double mean_local = m_Mean;

      if( m_bUseLocalThreshold ){
         if( rms_for_dm.size() > 0 && y < rms_for_dm.size() ){
            double rms = rms_for_dm[y];
            if( rms > 0 && !isnan(rms) ){
               printf("for Y = %d using rms = %.4f (instead of generic = %.4f) and threshold = %.4f ( generic = %.4f )\n",y,rms,m_RMS,(threshold_in_sigma*rms + m_Mean),threshold_jy);        
               threshold_jy = threshold_in_sigma*rms + m_Mean;         
               rms_local = rms;
            }
         }
      }
         
      for(int x=minX;x<maxX;x++){      
          if( m_bUseLocalThreshold ){
             if( p_rms_per_count && p_rms_per_count->size() > 0 && pCountMap ){
                // if RMS per count is available use with higher priority than RMS(DM) :
                int map_count = pCountMap->getXY( x, y );
                if( map_count < p_rms_per_count->size() ){
                   double rms = (*p_rms_per_count)[ map_count ];

                   // fitted : 2.137 / sqrt(map_count ) :
                   if ( false ){ // use fitted dependence ?
                      rms = 2.137 / sqrt(map_count );
                   }

                   rms_local = rms;                
                   threshold_jy = threshold_in_sigma*rms + m_Mean;
                }else{
                   printf("ERROR in code : pixel (%d,%d) was calculated out of %d pixels which is > rms_per_count.size() = %d\n",x,y,map_count,int(p_rms_per_count->size()));
                }
             }
          }
      
          double val = getXY(x,y);
          
          if( val > threshold_jy ){
             CSource* pCloseSource = m_Sources.find( x, y, 3 );
             if( x > debug_maxX ){
                debug_maxX = x;
             }
             
             if( pCloseSource ){
                if( val > pCloseSource->flux ){
                   pCloseSource->Update( x, y, val );
                }
             }else{
                m_Sources.Add( x, y, val, rms_local, mean_local );
             }             
          }
      }
   }

   printf("DEBUG : identified %d sources above the threshold (local = %d) , debug_maxX = %d\n",int(m_Sources.size()),m_bUseLocalThreshold,debug_maxX);

   string szOutTxtFile,szOutCandFile;
   change_ext( szOutRegFile, "txt", szOutTxtFile );   
   change_ext( szOutRegFile, "cand", szOutCandFile );
   MyOFile regfile( szOutRegFile ,"w");
   MyOFile txtfile( szOutTxtFile.c_str(), "w" );
   MyOFile candfile( szOutCandFile.c_str(), "w");
   printf("Sources in reference image %s are :\n",GetFileName());

   // # S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd
   // 10.01 604 0.6040 31 176 60.03 0 60096.430354213   
   candfile.Printf("# S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, DTM, X, Y\n");
   
   if( bFluxInSigmas ){
      txtfile.Printf("# X   Y    SNR N_pix\n");
   }else{
      txtfile.Printf("# X   Y    FLUX[Jy] N_pix\n");
   }
   for(int s=0;s<m_Sources.size();s++){
      CSource& source = m_Sources[s];
      int n_pix = 0;
      if( pCountMap ){
         n_pix = pCountMap->getXY( (int)source.x ,(int)source.y );
      }

   
      double snr = (source.flux - source.mean) / source.rms; // m_Mean added
      if( bFluxInSigmas ){
         printf("\t(%.1f,%.1f) %.2f snr\n",source.x,source.y,snr);
         regfile.Printf("\tcircle %.1f %.1f 3 # color=green snr = %.5f , n_pix = %d , snr = (%.4f - %.4f) / %.4f\n",source.x,source.y,snr,n_pix,source.flux,source.mean,source.rms);
         txtfile.Printf("%.1f %.1f %.4f %d\n",source.x,source.y,snr,n_pix);
      }else{
         printf("\t(%.1f,%.1f) %.2f Jy\n",source.x,source.y,source.flux);
         regfile.Printf("\tcircle %.1f %.1f 3 # color=green flux = %.5f Jy , n_pix = %d\n",source.x,source.y,source.flux,n_pix);
         txtfile.Printf("%.1f %.1f %.4f %d\n",source.x,source.y,source.flux,n_pix);         
      }
   
      // # S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd
      // 10.01 604 0.6040 31 176 60.03 0 60096.430354213   
      int x=-1,y=-1;
      sscanf(GetFileName(),"%d_%d_series.fits",&x,&y);
      candfile.Printf("%.2f %.2f %.6f - - %.3f %s 0.00 - %d %d\n",snr,source.x,source.x*cdelt1,source.y,GetFileName(),x,y);
   }   
   
   printf("DEBUG( CMWAFits::FindSources ) : found candidates %d above global threshold %.4f\n",int(m_Sources.size()),threshold_jy_global);
   
   return m_Sources.size();
}

int CMWAFits::FindSourcesSNR( double threshold_snr /* =5.00 */, int border /* = 20 */ , const char* szOutRegFile /* ="sources_snr.reg" */ , int minX /* =0 */, int maxX /* =-1 */ , CBgFits* pCountMap /*=NULL*/ , bool bSavePhysical /*=false*/ )
{
   printf("PROGRESS CMWAFits::FindSourcesSNR size %d x %d, PARAMS : border = %d\n",int(GetXSize()),int(GetYSize()),border);

   m_Sources.clear();

   if( maxX <= 0 ){
      maxX = m_SizeX;
   }
   
   for(int y=0;y<m_SizeY;y++){
      for(int x=minX;x<maxX;x++){
          double val = getXY(x,y);
          
          if( val > threshold_snr ){
             CSource* pCloseSource = m_Sources.find( x, y, 3 );
             
             if( pCloseSource ){
                if( val > pCloseSource->flux ){
                   pCloseSource->Update( x, y, val );
                }
             }else{
                m_Sources.Add( x, y, val );
             }             
          }
      }
   }

   string szOutTxtFile,szOutTxtFile_DTS;
   change_ext( szOutRegFile, "txt", szOutTxtFile );   
   change_ext( szOutRegFile, "_dts.txt", szOutTxtFile_DTS, true );   
   MyOFile regfile( szOutRegFile ,"w");
   MyOFile txtfile( szOutTxtFile.c_str(), "w" );
   MyOFile txtfile_dts( szOutTxtFile_DTS.c_str(), "w" );
   printf("Sources in reference image %s are :\n",GetFileName());

   if( bSavePhysical ){    
      txtfile.Printf("# SAMPNO_DYNASPEC  DM  SNR  N_PIX\n");
      txtfile_dts.Printf("# SAMPNO_DTS  DM  SNR  N_PIX\n");
   }else{
      txtfile.Printf("# X   Y    SNR  N_PIX\n");      
   }
   for(int s=0;s<m_Sources.size();s++){
      CSource& source = m_Sources[s];
      
      int n_pix = 0;
      if( pCountMap ){
         n_pix = pCountMap->getXY( (int)source.x ,(int)source.y );
      }

      printf("\t(%.1f,%.1f) snr = %.2f n_pix = %d (FindSourcesSNR)\n",source.x,source.y,source.flux,n_pix);      
   
      if( bSavePhysical ){   
         ReadCRValues();
         
         double x_value = get_x_physical_value(source.x);
         double y_value = get_y_physical_value(source.y);
         
         if( x_value >= border && x_value <= (GetXSize()-border) ){ // ignore negative arrival times (just noise)            
            double sampno = source.x - fabs(crval1/cdelt1);
            txtfile.Printf("%.1f %.1f %.4f %d %.1f %.1f\n",sampno,y_value,source.flux,n_pix,source.x,source.y);
            txtfile_dts.Printf("%.1f %.1f %.4f %d %.1f %.1f\n",source.x,y_value,source.flux,n_pix,source.x,source.y);
            regfile.Printf("\tcircle %.1f %.1f 3 # color=green snr = %.5f , n_pix = %d , time = %.4f [sec] (FindSourcesSNR)\n",source.x,source.y,source.flux,n_pix,x_value);
         }else{
            printf("WARNING : candidate at (x,y) = (%.1f,%.1f) skipped due to physical values (%.4f,%.4f) not within allowed range\n",source.x,source.y,x_value,y_value);
         }
      }else{         
         txtfile.Printf("%.1f %.1f %.4f %d\n",source.x,source.y,source.flux,n_pix);
         regfile.Printf("\tcircle %.1f %.1f 3 # color=green snr = %.5f , n_pix = %d (FindSourcesSNR)\n",source.x,source.y,source.flux,n_pix);
      }
   }
   
   return m_Sources.size();

}

int CMWAFits::CopyWithoutBadColumns( CMWAFits& out_fits, vector<int>& good_columns )
{
   out_fits.Realloc( good_columns.size(), GetYSize(), FALSE );
   for(int i=0;i<good_columns.size();i++){
       double column = good_columns[i];
       
       for(int y=0;y<GetYSize();y++){
          out_fits.setXY( column, y, getXY(column,y) );
       }
   }

   return good_columns.size();
}

int CMWAFits::CutValues( double threshold_in_sigma, vector<int>& good_columns, double zero_value, int rms_border )
{
   good_columns.clear();
   printf("COLUMNS :\n");
   for(int x=0;x<GetXSize();x++){
       double rms_column = 0.00;
       double mean_column = MeanColumn( x , &rms_column );
       
       printf("\t%d : mean = %.5f , rms = %.5f\n",x,mean_column,rms_column);
       
       if( fabs(mean_column) < 0.000000001 ){
           for(int y=0;y<GetYSize();y++){
               setXY(x,y,zero_value);
           }
           printf("\tColumn %d is BAD !!!\n",x);
       }else{
           good_columns.push_back( x );
           printf("\tColumn %d is good\n",x);
       }
   }
   
   double sum=0.00,sum2=0.00;
   int count=0;
   for(int y=rms_border;y<(GetYSize()-rms_border);y++){
       for(int x=rms_border;x<(GetXSize()-rms_border);x++){
           double val = getXY(x,y);
           
           if( val > -100 ){
               sum += val;
               sum2 += val*val;
               count++;
           }
       }
   }
   
   m_Mean = sum / count;
   m_RMS  = sqrt(sum2/count - m_Mean*m_Mean);
   
   double threshold = m_Mean + threshold_in_sigma*m_RMS;
   
   printf("CutValues : mean = %.5f , rms = %.5f , threshold = %.5f\n",m_Mean,m_RMS,threshold);
   
   count=0;
   for(int y=0;y<GetYSize();y++){
       for(int x=0;x<GetXSize();x++){
           double val = getXY(x,y);
           
           if( val < threshold )
           {
              setXY(x,y,zero_value);           
           }else{
              // set values in SIGMAS :
              // setXY(x,y, (val/m_RMS) );
              setXY( x,y, val );
              count++;
           }
       }
   }
   
   printf("CutValues : accpeted %d our of %d pixels\n",count,int(GetXSize()*GetYSize()));
   
   return count;
}


CMWADataCube::CMWADataCube()
: m_Obsid(-1), m_Timesteps(0), m_Channels(0), m_Metafits(NULL), m_pMeanImage(NULL), m_pMedianImage(NULL), m_pRMSImage(NULL), m_FreqCenterMHz(-1)
{

}

CMWADataCube::CMWADataCube( int obsid, int n_channels, int timesteps, int start_timeindex )
: m_Obsid(obsid), m_Timesteps(timesteps), m_Channels(n_channels), m_Metafits(NULL), m_pMeanImage(NULL), m_pMedianImage(NULL), m_pRMSImage(NULL),
  m_FirstCorrectTimestep(-1), m_FirstCorrectChannel(-1), m_StartTimeIndex(start_timeindex), m_FirstCoarseChannel(133)
{
   printf("DEBUG : contructing CMWADataCube::CMWADataCube object\n");fflush(stdout);   

   // 2021-05-03 - new implementation 
   Init( obsid, n_channels, timesteps, start_timeindex );   
}

void CMWADataCube::Init( int obsid, int n_channels, int timesteps, int start_timeindex, double freq_center_mhz, double freq_lower_mhz, double freq_upper_mhz ) 
{
   m_Obsid     = obsid;
   m_Timesteps = timesteps;
   m_Channels  = n_channels;
   m_StartTimeIndex = start_timeindex;
   m_FreqLowerMHz = freq_lower_mhz;
   m_FreqUpperMHz = freq_upper_mhz;
   m_FreqCenterMHz = freq_center_mhz;
   
   if( m_Timesteps <= 0 ){
      // automatically check number of timesteps :
      printf("AUTO-ACTION : need to automatically calculate number of timesteps\n");fflush(stdout);
      m_Timesteps = count_time_steps();      
   }

   if( m_FreqCenterMHz <= 0 ){
      m_FreqCenterMHz = 156*MWA_COARSE_CHANNEL;
   }   
   if( m_FreqLowerMHz <= 0 ){
      m_FreqLowerMHz = m_FreqCenterMHz - MWA_COARSE_CHANNEL/2.00; // default for center channel 145 ( 133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156 )
   }
   if( m_FreqUpperMHz <= 0 ){
      m_FreqUpperMHz = m_FreqCenterMHz + MWA_COARSE_CHANNEL/2.00; // default for center channel 145 ( 133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156 )
   }
   
   // init with default coarse channels as above :
   int default_coarse_channels[] = { 133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156 };
   for(int ch=0;ch<24;ch++){
       m_CoarseChannels.push_back( default_coarse_channels[ch] );
   }

   printf("Creating timesteps vector of size = %d\n",m_Timesteps);fflush(stdout);
   vector< CMWAFits* > time_steps_template;
   time_steps_template.assign( m_Timesteps, NULL );

   assign( n_channels, time_steps_template );
   
   printf("Allocated data cube of CHANNELS x TIMESTEPS = %d x %d\n",((int)this->size()),((int)(*this)[0].size()));fflush(stdout);
   
   m_StartUnixTime = gps2ux( m_Obsid );

   // init flags structure :
   vector< int > flags_template;
   flags_template.assign( m_Timesteps, 0 );
   m_ImageFlags.assign( n_channels, flags_template );

   DebugDump();
   
}

CMWADataCube::~CMWADataCube()
{
   if( m_Metafits ){
       delete m_Metafits;
   }
   
   if( m_pMeanImage ){
       delete m_pMeanImage;
   }
   if( m_pRMSImage ){
       delete m_pRMSImage;
   }
   if( m_pMedianImage ){
       delete m_pMedianImage;
   }

   for(int t=0;t<m_Timesteps;t++){
       for(int ch=0;ch<m_Channels;ch++){
           if( ((*this)[ch][t]) ){
               delete ((*this)[ch][t]);
           }
       }
   }
}

int CMWADataCube::GetDynamicSpectrum( int x, int y, CBgFits& out_dynaspec )
{
    int xSize = 0, ySize = 0;

    // set ZEROs to avoid "crazy" values     
    out_dynaspec.SetValue( 0 );

    if( this->size()>0 && ((*this)[0]).size()>0 ){
       CMWAFits* pImage = GetImage( m_FirstCorrectChannel , m_FirstCorrectTimestep );
       if( pImage ){
//           xSize = ((*this)[0][0])->GetXSize();
//           ySize = ((*this)[0][0])->GetYSize();
           xSize = pImage->GetXSize();
           ySize = pImage->GetYSize();
       }else{
           printf("ERROR in code CMWADataCube::GetDynamicSpectrum pImage = NULL\n");
           exit(-1);           
       }
    }
        
    if( x>=0 && y>=0 && x<xSize && y<ySize ){
       if( out_dynaspec.GetXSize() != m_Timesteps || out_dynaspec.GetYSize() != m_Channels ){
          out_dynaspec.Realloc( m_Timesteps, m_Channels, FALSE );
       }
       for(int t=0;t<m_Timesteps;t++){
           for(int ch=0;ch<m_Channels;ch++){
               if( ((*this)[ch][t]) ){
                   double val = ((*this)[ch][t])->getXY( x, y );
                   out_dynaspec.setXY( t, ch, val );
                   if( gDebugLocal > 1 ){
                      printf("DEBUG : (%d,%d) = %.8f\n",t,ch,val);                                     
                   }
               }else{
                   printf("WARNING : no or flagged data for channel = %d , timestep = %d\n",ch,t);
               }
           }
       }
    }else{
       printf("ERROR : pixel (%d,%d) is outside the image of size %d x %d pixels\n",x,y,xSize,ySize);
    }

    return -1;
}

int CMWADataCube::GetOptimalImage( double unixtime, double& min_dist )
{   
   int closest_index = GetClosestImage( unixtime, min_dist );
   if (  closest_index > 0 && closest_index < (*this)[0].size() ){
      double unixtime_closest  = ((*this)[0][closest_index])->GetUnixTime();   
      if( ((*this)[0][closest_index-1]) ){
         double unixtime_previous = ((*this)[0][closest_index-1])->GetUnixTime();      
         double middle = (unixtime_closest + unixtime_previous)/2.00;
      
         if( unixtime >= middle && unixtime < unixtime_closest && unixtime > unixtime_previous ){
            // return previous index to include more FRB flux :
            return (closest_index-1);
         }
      }
   }
   
   return closest_index;
}

int CMWADataCube::GetClosestImage( double unixtime, double& min_dist )
{
   min_dist = 1e20;
   int out_index   = -1;
   
   for(int t=0;t<m_Timesteps;t++){
      if( this->size() > 0 && t < (*this)[0].size() ){
         if( ((*this)[0][t]) ){
            double fits_unixtime = ((*this)[0][t])->GetUnixTime();
            double dist = fabs(fits_unixtime-unixtime);
            if( gDebugLocal > 2 ){
               HeaderRecord* pKeyword = ((*this)[0][t])->GetKeyword("DATE-OBS");
               string szDate="???";
               if( pKeyword ){
                  szDate = pKeyword->Value.c_str();
               }
               printf("\tDEBIG : CMWADataCube::GetClosestImage(%.2f) %d = %.2f( %s ) -> dist = %.4f sec ( min_dist = %.2f sec )\n",unixtime,t,fits_unixtime,szDate.c_str(),dist,min_dist);
            }
      
            if( dist < min_dist ){
               min_dist = fabs(fits_unixtime-unixtime);
               out_index = t;
            }
         }else{
            if( gDebugLocal>2 ){ printf("WARNING : Null pointed at timeindex = %d , channel = 0 -> ignored\n",t); }
         }
      }else{
         printf("ERROR : m_Timestamps = %d, this->size() = %d, (*this)[0].size() = %d < t = %d\n",m_Timesteps,int(this->size()),int((*this)[0].size()),t);
      }
   }
   
   if( out_index >= 0 && out_index < m_Timesteps ){
//      return ((*this)[0][out_index])->GetUnixTime();
      if( gDebugLocal>=2 ){printf("\tRETURN %d\n",out_index);}
      return out_index;
   }
   
   return -1;
}

CMWAFits* CMWADataCube::GetImage( int channel, int timeindex )
{
    return ((*this)[channel][timeindex]);
}

int CMWADataCube::GetXSize() // { return GetImage(0,0)->GetXSize(); }
{
   CMWAFits* pImage = GetImage(0,0);
   if( pImage ){
      return pImage->GetXSize();
   }   
   printf("ERROR : no image in the data cube !\n");
   
   return -1;
}

int CMWADataCube::GetYSize() // { return GetImage(0,0)->GetYSize(); }
{
   CMWAFits* pImage = GetImage(0,0);
   if( pImage ){
      return pImage->GetYSize();
   }   
   printf("ERROR : no image in the data cube !\n");
   
   return -1;
}



void CMWADataCube::DebugDump( eDumpType dump_type /* = eDumpType_Ptr */ )
{
   printf("CMWADataCube::DebugDump : gDebug = %d\n",gDebugLocal);
   if( gDebugLocal>0 ){   
       printf("Test dump of data cube\n");
       for(int t=0;t<m_Timesteps;t++){
          printf("Timeindex%03d : ",t);
          for(int ch=0;ch<m_Channels;ch++){
             switch( dump_type ){
                case eDumpType_Ptr :
                    printf("%p ",(void *)((*this)[ch][t]));         
                    break;
                
                case eDumpType_CenterPixel :
                    if( ((*this)[ch][t]) ){
                        printf("%.3f ",((*this)[ch][t])->valXY( ((*this)[ch][t])->GetXSize()/2, ((*this)[ch][t])->GetYSize()/2 ) );
                    }else{
                        printf("%p ",((*this)[ch][t]));
                    }
                    break;
                    
                default :
                    printf("%p ",((*this)[ch][t]));
                    break;
             }
          }
          printf("\n");      
       }
   }   
}

int  CMWADataCube::DoesDirExist( const char* szDirPath )
{
    int bExists = 0;
    struct stat s;
    int err = stat( szDirPath, &s);
    if(-1 == err) {
        if(ENOENT == errno) {
            /* does not exist */
        } else {
            perror("stat");
            exit(1);
        }
    } else {
        if(S_ISDIR(s.st_mode)) {
            /* it's a dir */
            bExists = 1;
        } else {
            /* exists but is no dir */
       }
    }
    
    return bExists;
}    
    

int CMWADataCube::count_time_steps( const char* dir_template, int max_timesteps /*=1000*/ )
{
    char szDir[1024];
    int timesteps = 0;
    int max_existing_timestep = -1;
    
    while( timesteps < max_timesteps ){
       sprintf(szDir,dir_template,timesteps);
       
       if( DoesDirExist(  szDir ) > 0 ){
          max_existing_timestep = timesteps;
       }
       timesteps++;
    }
    printf("AUTO-ACTION : determined number of time steps to be %d\n",max_existing_timestep);
 
    return max_existing_timestep;       
}


int CMWADataCube::generate_filfile( int x, int y, const char* out_filename, double f_start, double foff, double uxstart, double bw, double inttime, bool bSaveRescaledFits )
{
#ifdef _USE_FILFILE_
   printf("INFO : generating filfile at pixel (%d,%d) fch1 = %.4f MHz, foff = %.4f MHz, uxstart = %.8f, bw = %.4f MHz , save_rescaled_fits = %d\n",x,y,f_start,foff,uxstart,bw,bSaveRescaledFits);

   // const char* change_ext(const char* name,const char* new_ext,string& out, bool bChangeFileName=false);
   string szOutFitsFile, szOutSpecFile;
   char szOutFileName[1024];
   if( out_filename ){
      strcpy( szOutFileName , out_filename );                  
   }else{
      sprintf(szOutFileName,"filterbank_atpix_%05d_%05d.fil",x,y);
      out_filename = szOutFileName;
   }
/*   if( out_specfilename ){
      strcpy( szOutSpecFile , out_specfilename );                  
   }else{
      sprintf(szOutSpecFile,"filterbank_atpix_%05d_%05d_spectrum.txt",x,y);
      out_filename = szOutSpecFile;
   }*/
   change_ext( szOutFileName, "fits", szOutFitsFile, 0 ); 
   change_ext( szOutFileName, "_spectrum.txt", szOutSpecFile, 1);
   
   if( x < 0 ){
      x = m_Metafits->GetXSize() / 2;
   }
   if( y < 0 ){
      y = m_Metafits->GetYSize() / 2;
   }
   
   CBgFits out_dynaspec( m_Timesteps, m_Channels );
   GetDynamicSpectrum( x, y, out_dynaspec );
   out_dynaspec.start_freq = f_start;
   out_dynaspec.delta_freq = foff;
   out_dynaspec.dtime_fs   = int(uxstart);  
   out_dynaspec.dtime_fu   = (uxstart - int(uxstart))*1000000.00;
   // void PrepareBigHornsHeaderTransposed( double ux_start, double _inttime, double freq_start, double delta_freq_mhz )
   double center_of_lowest_channel_mhz = f_start-bw+fabs(foff)/2.00; // to put into CRVAL1 or 2 
   out_dynaspec.PrepareBigHornsHeaderTransposed( uxstart, inttime, center_of_lowest_channel_mhz,  fabs(foff) );

   // calculate min and max for re-scaling for fil file:
   // double GetStat( double& mean, double& rms, double& minval, double& maxval, int x_start=0, int y_start=0, int x_end=-1, int y_end=-1 );
   double mean,rms,minval,maxval;
   out_dynaspec.GetStat( mean, rms, minval, maxval );
   
   bool b_rescale_inside_fits2fil = true;
   if( bSaveRescaledFits ){
      // FITS dynamic spectrum is re-scaled :
      int xSize = out_dynaspec.GetXSize();
      int ySize = out_dynaspec.GetYSize();
      for(int x=0;x<xSize;x++){
         for(int y=0;y<ySize;y++){
            double val = out_dynaspec.getXY(x,y);
            float scaled_value = val;
         
            scaled_value = 255.00*(scaled_value - minval)/(maxval-minval);

            if( fabs(maxval) > 1.00 ){
              if( scaled_value > maxval ){
                 scaled_value = maxval;
               }
            }
            out_dynaspec.setXY(x,y,scaled_value);
         }
      }
      b_rescale_inside_fits2fil = false;
      
      // write for debuging :   
      out_dynaspec.WriteFits( szOutFitsFile.c_str() );
   }else{   
      // write for debuging :   
      out_dynaspec.WriteFits( szOutFitsFile.c_str() );
   
      // rescaling inside fits2fil (otherwise if not set -> no rescaling inside fits2fil :
      CFilFile::gMinFITS_Value = minval;
      CFilFile::gMaxFITS_Value = maxval;
      b_rescale_inside_fits2fil = true;
   }
   printf("CMWADataCube::generate_filfile GetStat on dynamic spectrum found values range %.8f - %.8f\n",minval,maxval);
   
   
   
   
   
   // fits2fil( CBgFits& infits, const char* szOutFilFile, const char* szOutSpecFile, int n_channels, int n_timesteps );
   CFilFile::fits2fil( out_dynaspec, szOutFileName, szOutSpecFile.c_str(), m_Channels, m_Timesteps , b_rescale_inside_fits2fil );
   
   // SigprocFile( int nbits, int nifs, int nchans, double fch1, double foff, double tstart, double tsamp );
//   SigprocFile filfile( 8, 1, m_Channels, 
//   CFilFile out_filfile( szOutFileName );
//   out_filfile.ParseHeader( *m_Metafits );
#else
   printf("ERROR : function CMWADataCube::generate_filfile is not implemented !\n");
   return -1;
#endif

   return 1;
}


bool CMWADataCube::Avg( CMWADataCube& right )
{
   return Operation( right, eAvgImages );
}

bool CMWADataCube::Operation( CMWADataCube& right , eCalcFitsAction_T oper )
{
   if( oper != eAvgImages ){
      printf("ERROR : unknown operation %d\n",oper);
      return false;
   }

   // first check dimenssions - must be the same :
   if( m_Channels != right.m_Channels ){
      printf("ERROR(CMWADataCube::Avg) : number of channels are not equal ( %d != %d )\n",m_Channels,right.m_Channels);
      return false;
   }
   if( m_Timesteps != right.m_Timesteps ){
      printf("ERROR(CMWADataCube::Avg) : number of timesteps are not equal ( %d != %d )\n",m_Timesteps,right.m_Timesteps);
      return false;
   }
   
   // compare files names (after removing XX and YY) - should be the same
   for(int t=0;t<m_Timesteps;t++){
      for(int ch=0;ch<m_Channels;ch++){
         CMWAFits* pImageLeft = GetImage( ch, t );
         CMWAFits* pImageRight = right.GetImage( ch, t );
         
         if( !pImageLeft ){
            printf("ERROR(CMWADataCube::Avg) : no LEFT image in channel = %d and timestep = %d\n",ch,t);
            return false;
         }
         if( !pImageRight ){
            printf("ERROR(CMWADataCube::Avg) : no RIGHT image in channel = %d and timestep = %d\n",ch,t);
            return false;
         }
         
         
         const char* pFileNameLeft = pImageLeft->GetFileName();
         const char* pFileNameRight = pImageRight->GetFileName();
         
         if( pFileNameLeft && pFileNameRight ){
            if( strlen(pFileNameLeft) && strlen(pFileNameRight) ){
               char szFileNameLeft[1024],szFileNameRight[1024];
               
               strcpy(szFileNameLeft,pFileNameLeft);
               strcpy(szFileNameRight,pFileNameRight);
               
               char* ptrXX = (char*)strstr(szFileNameLeft,"XX");
               char* ptrYY = (char*)strstr(szFileNameRight,"YY");
               
               if( !ptrXX || !ptrYY ){
                  printf("WARNING : one of the files %s or %s does not have XX or YY in the file name !\n",szFileNameLeft,szFileNameRight);
               }
               
               // making XX and YY file names the same to use strcmp for comparison
               ptrXX[0] = '_';
               ptrXX[1] = '_';
               ptrYY[0] = '_';
               ptrYY[1] = '_';
               if( strcmp( szFileNameLeft,szFileNameRight) ){
                  printf("ERROR : different file names for channel = %d , timestep = %d |%s| != |%s|\n",ch,t,szFileNameLeft,szFileNameRight);
                  return false;
               }
            }
         }
         
         if( pImageLeft->GetXSize() != pImageRight->GetXSize() || pImageLeft->GetYSize() != pImageRight->GetYSize() ){
            printf("ERROR : image dimenssion not the same %d != %d or %d != %d at channel = %d , timestep = %d\n",int(pImageLeft->GetXSize()),int(pImageRight->GetXSize()),int(pImageLeft->GetYSize()),int(pImageRight->GetYSize()),ch,t);
            return false;
         }
         
         
         // eAvgImages
         
         if( oper == eAvgImages ){
            // everything is the same for both images, can do a loop overpixels and perform operation :
            int xSize = pImageLeft->GetXSize();
            int ySize = pImageLeft->GetYSize();
            for(int y=0;y<ySize;y++){
               for(int x=0;x<xSize;x++){
                  double val_left = pImageLeft->getXY(x,y);
                  double val_right = pImageRight->getXY(x,y);
                  double val_out = (val_left + val_right) / 2.00;
                  
                  pImageLeft->setXY( x, y, val_out );
               }
            }
         }
      }
   }
   
   // Stokes I = (x+y)/2
   // see test.cpp
   
   return true;
}

int CMWADataCube::SaveChannelFitsFiles( int n_images )
{
    if( n_images < 0 ){
       n_images = m_Timesteps;
    }
    
    if( this->size()>0 && ((*this)[0]).size()>0 ){
       CMWAFits* pImage = GetImage( m_FirstCorrectChannel , m_FirstCorrectTimestep );
       if( !pImage ){
           printf("ERROR in code CMWADataCube::GetDynamicSpectrum pImage = NULL\n");
           exit(-1);
       }
    }
        
    for(int t=0;t<n_images;t++){
        string szFirstFits="tmp.fits"; // -> out directory tmp/
        if( (*this)[0][t] && ((*this)[0][t])->GetFileName() ){
           szFirstFits = ((*this)[0][t])->GetFileName();
        }
        string szOutDir = "tmp";
               
        for(int ch=0;ch<m_Channels;ch++){
            if( ((*this)[ch][t]) ){
               // const char* change_ext(const char* name,const char* new_ext,string& out, bool bChangeFileName=false);
               change_ext( szFirstFits.c_str(), "", szOutDir, true ); 
//               printf("szOutDir = %s\n",szOutDir.c_str());
            }else{
                printf("WARNING : no or flagged data for channel = %d , timestep = %d\n",ch,t);
            }
               
            char szOutFits[1024];
            sprintf(szOutFits,"%s/%s_ch%03d.fits",szOutDir.c_str(),szOutDir.c_str(),ch);
               
            printf("Saved fits file : %s\n",szOutFits);
            ((*this)[ch][t])->SetKeys( m_Metafits-> GetKeys() );
            if( ((*this)[ch][t])->WriteFits( szOutFits ) ){
               printf("ERROR : while writting FITS file %s\n",szOutFits);                  
            }
        }
    }

    return 1;

}

int CMWADataCube::ReadFitsCube( const char* fits_list, double freq_lower_mhz, double freq_upper_mhz, int bAutoDetect, int bReadImage, int bIgnoreHeaderErrors )
{
   // read list of FITS files with data cubes :
   vector<string> vec_fits_list;
   if( bg_read_list( fits_list, vec_fits_list ) <= 0 ){
      printf("ERROR : could not read list of fits files with data cubes from file %s\n",fits_list);
      exit(-1);
   } 
   
   // initialise things normally done in the contructor :
   fitsfile *fp;
   int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
   int anaxis;
   long npixels = 1;

   
   // read first fits without data as METAFITS :
   m_Metafits = new CBgFits( vec_fits_list[0].c_str() );
   if( m_Metafits->ReadFits( vec_fits_list[0].c_str() , 1 , 0 ) ){
       printf("ERROR : could not read header of FITS file %s - this error can not be tolerated -> exiting now !!!\n",vec_fits_list[0].c_str());
       exit(-1);
   }
   int xSize  = m_Metafits->GetXSize();
   int ySize  = m_Metafits->GetYSize();
   int sizeXY = m_Metafits->GetXSize()*m_Metafits->GetYSize();   

   const char* fits_file = vec_fits_list[0].c_str();
   printf("INFO : reading FITS cube %s\n",fits_file);
   fits_open_file( &fp, fits_file, READONLY, &status ); /* open input images */
   if (status) {
      fits_report_error(stderr, status); /* print error message */
      return(status);
   }

   fits_get_img_dim(fp, &anaxis, &status);  /* read dimensions */
   printf("INFO : image  %s has %d axis\n",fits_file,anaxis);
   
   long* firstpixel = new long[anaxis];
   int image_type = TFLOAT;
   firstpixel[0] = 1; // fits.GetXSize()
   firstpixel[1] = 1; // fits.GetYSize()
   firstpixel[3] = 1;
   

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
         sprintf(szTmp,"%d x ",int(anaxes[i]));
      }else{
         sprintf(szTmp,"%d",int(anaxes[i]));
      }
      szDim += szTmp;
   }
   printf("INFO : image dimensions = %s\n",szDim.c_str());

   double gpsTime = m_Metafits->dtime_fs + m_Metafits->dtime_fu / 1000000.00 - 315964783; // convert UX time to GPS time   
   int n_channels = anaxes[2];
   int n_timesteps = vec_fits_list.size();
   m_Channels = n_channels;
   m_Timesteps = n_timesteps;
   double freq_center_mhz = ( freq_lower_mhz + freq_upper_mhz ) / 2.00;
   Init( gpsTime, n_channels, n_timesteps, 0, freq_center_mhz, freq_lower_mhz, freq_upper_mhz );
   
   fits_close_file(fp, &status);
   if (status) {
      fits_report_error(stderr, status); /* print error message */
      return(status);
   }
      
   for(int timeindex=0;timeindex<vec_fits_list.size();timeindex++){   
      fits_file = vec_fits_list[timeindex].c_str();
      
      // open FITS file 
      printf("PROGRESS : reading FITS file %s ...\n",fits_file);      
      fits_open_file( &fp, fits_file, READONLY, &status ); /* open input images */
      if (status) {
         fits_report_error(stderr, status); /* print error message */
         return(status);
      }

   
      // test reading all channels :
      // float* data = new float[sizeXY];
      
      for(int channel=0;channel<n_channels;channel++){
         firstpixel[2] = channel + 1; // in FORTRAN 1-indexed
         
         CMWAFits* pImage = GetImage( channel, timeindex );
         if( !pImage ){
            pImage = new CMWAFits( fits_file, xSize, ySize );            
            SetImage( channel, timeindex, pImage );
         }
         
         if( pImage ){
            pImage->SetFileName( fits_file );
         
            float* data = pImage->get_data();
            fits_read_pix( fp, image_type, firstpixel, sizeXY, NULL, data, NULL, &status); 
            if( gDebugLocal >= 1 ){
               printf("DEBUG1 : read image for channel = %d, timestamps = %d\n",channel,timeindex);
            }
         }else{
            printf("ERROR in code : could not find CBgFits image for channel = %d and timeindex = %d\n",channel, timeindex );
         }
      }

      // close 
      fits_close_file(fp, &status);
      if (status) {
         fits_report_error(stderr, status); /* print error message */
         return(status);
      }
   }
   m_FirstCorrectTimestep = 0;
   m_FirstCorrectChannel  = 0;
   printf("PROGRESS : finished reading FITS files from list file %s\n",fits_list);fflush(stdout);
  

   return 0;
}

// reading in channel (00000/ etc) and Time directory structures, it requires fits_list in each channel directory to exist;        
int CMWADataCube::ReadChanTime( const char* dir_template /*="%05d"*/,  // channel directory name
                                const char* image_template /*="dirty_image_%dT%d_real.fits"*/,
                                bool bExitOnReadError /*= true*/
                              )
{
   char szDir[1024],szFitsFilename[1024],szFullFitsPath[2048],szObsidFile[1024];
   int bAutoDetect=1;
   int x_size = -1;
   int y_size = -1;
   int read_images = 0;

   for(int ch=0;ch<m_Channels;ch++){
      std::vector<string> channel_fits_list;
      char szFitsList[128];
      sprintf(szFitsList,"%05d/fits_list",ch);
      
      if( bg_read_list( szFitsList , channel_fits_list ) > 0 ){   
         int n_timesteps = m_Timesteps;
         if( n_timesteps >= channel_fits_list.size() ){
            n_timesteps = channel_fits_list.size();
         }
      
         printf("DEBUG : channel %05d -> reading %d FITS files for specific timestamps\n",ch,n_timesteps);
         for(int t=0;t<n_timesteps;t++){ // m_Timesteps has to be known earlier see count_tim
            CMWAFits* pBgFits = new CMWAFits();
            if( !pBgFits ){
               printf("ERROR : could not allocate memory for object of class CMWAFits\n");
               exit(-1);
            }
            const char* szFitsName = channel_fits_list[ m_StartTimeIndex + t ].c_str();
   
            sprintf(szFitsFilename,"%s",szFitsName);
            sprintf(szFullFitsPath,"%05d/%s",ch,szFitsFilename);

            printf("Reading timestep = %d , channel = %d -> fits = %s , fits-obsid = %d\n",t,ch,szFullFitsPath,pBgFits->m_Obsid);

            if( !MyFile::DoesFileExist( szFullFitsPath ) ){
               mystring szFirstPath = szFullFitsPath;              
               mystring szGzippedFile = szFirstPath.c_str();
               szGzippedFile += ".gz";
               if ( MyFile::DoesFileExist( szGzippedFile.c_str() ) ){
                  printf("WARNING : reading compressed file %s\n",szGzippedFile.c_str());
                  sprintf(szFullFitsPath,"%s",szGzippedFile.c_str());                 
               }else{
                   sprintf(szFullFitsPath,"%s/wsclean_%d_timeindex%d-%04d-I-dirty.fits",szDir,m_Obsid,t,ch);                  
                   printf("WARNING : image %s could not be read -> trying old template %s ... \n",szFirstPath.c_str(),szFullFitsPath);
               }
            }

           
            if( MyFile::DoesFileExist( szFullFitsPath ) && (pBgFits->m_ReadStatus = pBgFits->ReadFits( szFullFitsPath, bAutoDetect ))==0 ){
                ((*this)[ch][t]) = pBgFits;
               
                x_size = pBgFits->GetXSize();
                y_size = pBgFits->GetYSize();               
               
                if( m_FirstCorrectTimestep < 0 && m_FirstCorrectChannel < 0 ){
                    m_FirstCorrectTimestep = t; // 20190822 - to only use m_StartTimeIndex for reading from the disk and writing results 
                    m_FirstCorrectChannel  = ch;
                    printf("AUTO-ACTION : first timestep with correctly read image is %d, channel = %d, image size = (%d,%d), this = %p\n",m_FirstCorrectTimestep,m_FirstCorrectChannel,x_size,y_size,this);
                }
                
                read_images++;
            }else{
                printf("\tERROR : file %s does not exist or could not be read\n",szFullFitsPath);
                delete pBgFits;
            }
         }
      }
   }

   if( read_images <= 0 ){
      if( bExitOnReadError ){
         printf("ERROR : could not read any image for obsID = %d -> exiting program\n",m_Obsid);
         exit(-1);      
      }else{
         printf("WARNING : could not read any image for obsID = %d -> returning from Read function (no exit required)\n",m_Obsid);
      }
   }
   
/*   if( m_Timesteps <= 0 ){
      printf("AUTO-ACTION : setting m_Timesteps := %d (was <=0 - which requires auto-action)\n",t);
      m_Timesteps = t;
   }*/
   
   DebugDump( eDumpType_CenterPixel );
   
/*   if( ! ReadMetaData( m_Obsid ) ){
      printf("ERROR : could not read metadata for obsid = %d\n",m_Obsid);
      exit(-1);
   }*/
   
//   printf("INFO : unixtime = %d\n",m_Metafits->dtime_fs); 
   printf("INFO : unixtime = %.2f (from gpstime = %d)\n",m_StartUnixTime,m_Obsid);

   if( x_size >0 && y_size>0 ){   
       m_pMeanImage = new CMWAFits( "mean.fits", x_size, y_size );
       m_pRMSImage = new CMWAFits( "rms.fits", x_size, y_size );
       m_pMedianImage = new CMWAFits( "median.fits", x_size, y_size );
   }
   
   return read_images;
}                              


int CMWADataCube::ReadBlinkImages( const char* image_template, /*="wsclean_%d_timeindex%03d-%04d-I-dirty.fits"*/
                     double time_resolution,
                     int n_seconds,
                     int first_second,
                     int first_coarse_channel, 
                     int n_coarse_channels,
                     const char* szDir, /* = "./" */
                     bool bExitOnReadError /* = true */ )
{
   char szFitsFilename[1024],szFullFitsPath[2048],szObsidFile[1064];
   int bAutoDetect=1;
   int x_size = -1;
   int y_size = -1;
   int read_images = 0;

   int time_steps_per_second = int(round(1.00/time_resolution));
   
   int newTimesteps = (n_seconds*time_steps_per_second);
   printf("DEBUG : ReadBlinkImages time_steps_per_second = %d -> total time steps = %d (was %d)\n",time_steps_per_second,newTimesteps,m_Timesteps);
   m_Timesteps = newTimesteps;

   int second=first_second; // TODO : needs to increase after full second added !!!
   for(int second=first_second;second<(first_second+n_seconds);second++){
      for(int t=0;t<(m_Timesteps);t++){ // 50 for 20ms 
         int total_fine_channel = 0;
         for(int cc=first_coarse_channel;cc<(first_coarse_channel+24);cc++){
             for(int ch=0;ch<m_Channels;ch++){
                sprintf(szObsidFile,"obsid.txt");
       
                 CMWAFits* pBgFits = new CMWAFits();
                 if( !pBgFits ){
                    printf("ERROR : could not allocate memory for object of class CMWAFits\n");
                    exit(-1);
                 }
           
                 pBgFits->m_Obsid = m_Obsid;
                 if( MyFile::DoesFileExist( szObsidFile )) {
                    MyFile obsid_file( szObsidFile );
                    const char* pLine = obsid_file.GetLine(TRUE);
                    pBgFits->m_Obsid = atol( pLine );              
                    m_Obsid = pBgFits->m_Obsid;
                    printf("DEBUG : read obsid = %d from fitsfile\n",pBgFits->m_Obsid);
                 } 
            
                 sprintf(szFitsFilename,image_template,second,t,cc,ch);
                 sprintf(szFullFitsPath,"%s/%s",szDir,szFitsFilename);
           
                 printf("Reading timestep = %d , coarse channel %d, channel = %d (total fine channel = %d)-> fits = %s , fits-obsid = %d\n",t,cc,ch,total_fine_channel,szFullFitsPath,pBgFits->m_Obsid);
           
                 if( !MyFile::DoesFileExist( szFullFitsPath ) ){
                    mystring szFirstPath = szFullFitsPath;              
                    mystring szGzippedFile = szFirstPath.c_str();
                    szGzippedFile += ".gz";
                    if ( MyFile::DoesFileExist( szGzippedFile.c_str() ) ){
                       printf("WARNING : reading compressed file %s\n",szGzippedFile.c_str());
                       sprintf(szFullFitsPath,"%s",szGzippedFile.c_str());                 
                    }else{
                       sprintf(szFullFitsPath,"%s/wsclean_%d_timeindex%d-%04d-I-dirty.fits",szDir,m_Obsid,t,ch);                  
                       printf("WARNING : image %s could not be read -> trying old template %s ... \n",szFirstPath.c_str(),szFullFitsPath);
                    }
                 }

           
                 if( MyFile::DoesFileExist( szFullFitsPath ) && (pBgFits->m_ReadStatus = pBgFits->ReadFits( szFullFitsPath, bAutoDetect ))==0 ){
                     ((*this)[ch][t]) = pBgFits;
               
                     x_size = pBgFits->GetXSize();
                     y_size = pBgFits->GetYSize();               
               
                     if( m_FirstCorrectTimestep < 0 && m_FirstCorrectChannel < 0 ){
                         m_FirstCorrectTimestep = t; // 20190822 - to only use m_StartTimeIndex for reading from the disk and writing results 
                         m_FirstCorrectChannel  = ch;
                         printf("AUTO-ACTION : first timestep with correctly read image is %d, channel = %d, image size = (%d,%d), this = %p\n",m_FirstCorrectTimestep,m_FirstCorrectChannel,x_size,y_size,this);
                     }
                 
                     read_images++;
                 }else{
                     printf("\tERROR : file %s does not exist or could not be read\n",szFullFitsPath);
                     delete pBgFits;
                 }
             }
          }
      }      
   }

   if( read_images <= 0 ){
      if( bExitOnReadError ){
         printf("ERROR : could not read any image for obsID = %d -> exiting program\n",m_Obsid);
         exit(-1);      
      }else{
         printf("WARNING : could not read any image for obsID = %d -> returning from Read function (no exit required)\n",m_Obsid);
      }
   }
   
/*   if( m_Timesteps <= 0 ){
      printf("AUTO-ACTION : setting m_Timesteps := %d (was <=0 - which requires auto-action)\n",t);
      m_Timesteps = t;
   }*/
   
   DebugDump( eDumpType_CenterPixel );
   
   if( !ReadMetaData( m_Obsid ) ){
      printf("ERROR : could not read metadata for obsid = %d\n",m_Obsid);
      exit(-1);
   }
   
//   printf("INFO : unixtime = %d\n",m_Metafits->dtime_fs); 
   printf("INFO : unixtime = %.2f (from gpstime = %d)\n",m_StartUnixTime,m_Obsid);

   if( x_size >0 && y_size>0 ){   
       m_pMeanImage = new CMWAFits( "mean.fits", x_size, y_size );
       m_pRMSImage = new CMWAFits( "rms.fits", x_size, y_size );
       m_pMedianImage = new CMWAFits( "median.fits", x_size, y_size );
   }
   
   return read_images;

}                     

int CMWADataCube::Read( const char* dir_template,  /* ="wsclean_timeindex%03d" */
                        const char* image_template, /* ="wsclean_%d_timeindex%03d-%04d-I-dirty.fits"*/ 
                        bool bExitOnReadError       /* = true */
                      )
{
   char szDir[1024],szFitsFilename[1024],szFullFitsPath[2048],szObsidFile[1064];
   int bAutoDetect=1;
   int x_size = -1;
   int y_size = -1;
   int read_images = 0;

   for(int t=0;t<(m_Timesteps);t++){ // m_Timesteps has to be known earlier see count_time_steps 
       sprintf(szDir,dir_template,t+m_StartTimeIndex);
   
       for(int ch=0;ch<m_Channels;ch++){
          sprintf(szObsidFile,"%s/obsid.txt",szDir);
       
           CMWAFits* pBgFits = new CMWAFits();
           if( !pBgFits ){
              printf("ERROR : could not allocate memory for object of class CMWAFits\n");
              exit(-1);
           }
           
           pBgFits->m_Obsid = m_Obsid;
           if( MyFile::DoesFileExist( szObsidFile )) {
              MyFile obsid_file( szObsidFile );
              const char* pLine = obsid_file.GetLine(TRUE);
              pBgFits->m_Obsid = atol( pLine );              
              m_Obsid = pBgFits->m_Obsid;
              printf("DEBUG : read obsid = %d from fitsfile\n",pBgFits->m_Obsid);
           } 
           
           sprintf(szFitsFilename,image_template,m_Obsid,t+m_StartTimeIndex,ch);
           sprintf(szFullFitsPath,"%s/%s",szDir,szFitsFilename);
           
           printf("Reading timestep = %d , channel = %d -> fits = %s , fits-obsid = %d\n",t,ch,szFullFitsPath,pBgFits->m_Obsid);
           
           if( !MyFile::DoesFileExist( szFullFitsPath ) ){
              mystring szFirstPath = szFullFitsPath;              
              mystring szGzippedFile = szFirstPath.c_str();
              szGzippedFile += ".gz";
              if ( MyFile::DoesFileExist( szGzippedFile.c_str() ) ){
                 printf("WARNING : reading compressed file %s\n",szGzippedFile.c_str());
                 sprintf(szFullFitsPath,"%s",szGzippedFile.c_str());                 
              }else{
                  sprintf(szFullFitsPath,"%s/wsclean_%d_timeindex%d-%04d-I-dirty.fits",szDir,m_Obsid,t,ch);                  
                  printf("WARNING : image %s could not be read -> trying old template %s ... \n",szFirstPath.c_str(),szFullFitsPath);
              }
           }

           
           if( MyFile::DoesFileExist( szFullFitsPath ) && (pBgFits->m_ReadStatus = pBgFits->ReadFits( szFullFitsPath, bAutoDetect ))==0 ){
               ((*this)[ch][t]) = pBgFits;
               
               x_size = pBgFits->GetXSize();
               y_size = pBgFits->GetYSize();               
               
               if( m_FirstCorrectTimestep < 0 && m_FirstCorrectChannel < 0 ){
                   m_FirstCorrectTimestep = t; // 20190822 - to only use m_StartTimeIndex for reading from the disk and writing results 
                   m_FirstCorrectChannel  = ch;
                   printf("AUTO-ACTION : first timestep with correctly read image is %d, channel = %d, image size = (%d,%d), this = %p\n",m_FirstCorrectTimestep,m_FirstCorrectChannel,x_size,y_size,this);
               }
               
               read_images++;
           }else{
               printf("\tERROR : file %s does not exist or could not be read\n",szFullFitsPath);
               delete pBgFits;
           }
       }
   }

   if( read_images <= 0 ){
      if( bExitOnReadError ){
         printf("ERROR : could not read any image for obsID = %d -> exiting program\n",m_Obsid);
         exit(-1);      
      }else{
         printf("WARNING : could not read any image for obsID = %d -> returning from Read function (no exit required)\n",m_Obsid);
      }
   }
   
/*   if( m_Timesteps <= 0 ){
      printf("AUTO-ACTION : setting m_Timesteps := %d (was <=0 - which requires auto-action)\n",t);
      m_Timesteps = t;
   }*/
   
   DebugDump( eDumpType_CenterPixel );
   
   if( ! ReadMetaData( m_Obsid ) ){
      printf("ERROR : could not read metadata for obsid = %d\n",m_Obsid);
      exit(-1);
   }
   
//   printf("INFO : unixtime = %d\n",m_Metafits->dtime_fs); 
   printf("INFO : unixtime = %.2f (from gpstime = %d)\n",m_StartUnixTime,m_Obsid);

   if( x_size >0 && y_size>0 ){   
       m_pMeanImage = new CMWAFits( "mean.fits", x_size, y_size );
       m_pRMSImage = new CMWAFits( "rms.fits", x_size, y_size );
       m_pMedianImage = new CMWAFits( "median.fits", x_size, y_size );
   }
   
   return read_images;
}                        

/*
int CMWADataCube::ReadNextStep()
{
   // delete [] in ((*this)[ch][t]) , where t=0
   // move everything by 1 back : 1->0, 2->1, 3->2, 4->3, .... N -> N-1
   // wczytac do N nowy obrazek !
   // cos trzeba po updatedowac jeszcze ???
   // gotowe ???
}
*/

bool CMWADataCube::GenerateMetaData( CBgFits& dyna_spec, CDedispSearch& dedisp, int obsid /* = -1 */, double start_freq /* = -1*/, double delta_freq /* = -1 */ )
{
   printf("################################## CMWADataCube::GenerateMetaData MODULE ##################################\n");   
   printf("WARNING : the value in brackets is approximate only - channels should be float precision to make it exact for any data (TODO)\n");
   m_Metafits = new CBgFits( "GENERATED" );
   m_Metafits->inttime = dyna_spec.inttime;   
   
   printf("DEBUG parameters : inttime = %.4f [sec] , start_freq = %.4f [MHz]\n",m_Metafits->inttime,dyna_spec.start_freq);
   
   if( start_freq < 0 ){
      start_freq= dyna_spec.start_freq;
   }
   if( delta_freq < 0 ){
      delta_freq= dyna_spec.delta_freq;
   }

   MWA_COARSE_CHANNELS = dedisp.m_MWADataCube.m_Channels;
   MWA_COARSE_CHANNEL  = ( delta_freq );
   printf("DEBUG N_FREQ_CHANNELS = %d , DELTA_FREQ = %.6f [MHz]\n",MWA_COARSE_CHANNELS,MWA_COARSE_CHANNEL);
   // dedisp.m_MWADataCube.m_FirstCoarseChannel = (dyna_spec.start_freq / dyna_spec.delta_freq);
   dedisp.m_MWADataCube.m_CoarseChannels.clear();
   for(int i=0;i<dedisp.m_MWADataCube.m_Channels;i++){
      double freq_mhz = start_freq + i*delta_freq; // assuming start_freq is the centre of the first channel it calculates cetres of channels


      double coarse_channel = ( freq_mhz / delta_freq );  // 2023-06-25 : int -> double to be able to have .5 of a channel which is the case here , round removed
      dedisp.m_MWADataCube.m_CoarseChannels.push_back( coarse_channel );
      dedisp.m_MWADataCube.m_fCoarseChannelFreqMHz.push_back( freq_mhz );
      
      int cc = CDedispSearch::find_value_double( m_CoarseChannels, coarse_channel );

      printf("   %d : %.6f MHz -> coarse channel = %.4f ( %.6f MHz - WARNING ONLY APPROXIMATE VALUE ), check cc = %d\n",i,freq_mhz,coarse_channel,(coarse_channel*delta_freq),cc);
   }
   dedisp.m_MWADataCube.m_FirstCoarseChannel = dedisp.m_MWADataCube.m_CoarseChannels[0];
   printf("The following metadata has been obtained from the dynamic spectrum file:\n");
   printf("MWA_COARSE_CHANNELS  = %d\n",MWA_COARSE_CHANNELS);
   printf("MWA_COARSE_CHANNEL   = %.4f MHz\n",MWA_COARSE_CHANNEL);
   printf("First coarse channel = %d\n",dedisp.m_MWADataCube.m_FirstCoarseChannel);
   
   return true;   
}


bool CMWADataCube::ReadMetaData( int obsid /* = -1 */ )
{
   if( obsid <= 0 ){
      obsid = m_Obsid;
   }

   // read metafits file 
   char szMetaFitsPath[1024];
   sprintf(szMetaFitsPath,"%d.metafits",obsid);
   if (!MyFile::DoesFileExist( szMetaFitsPath ) ){
      printf("WARNING : metafits file %s does not exist -> trying to get it ...\n",szMetaFitsPath);
      
      char szCmd[1024];      
      char szMetaFitsPathTmp[1024];
//      sprintf(szCmd,"wget http://mwa-metadata01.pawsey.org.au/metadata/fits/?obs_id=%d -O %d.metafits",obsid,obsid);
      
      sprintf(szMetaFitsPathTmp,"../../%d.metafits",obsid);
      if( MyFile::DoesFileExist(szMetaFitsPathTmp) ){
         printf("File %s exist -> copying ...\n",szMetaFitsPathTmp);
         sprintf(szCmd,"cp ../../%d.metafits .",obsid);
         printf("Executing command : %s\n",szCmd);
         int ret_val = system(szCmd);            
      }else{
         sprintf(szCmd,"wget http://ws.mwatelescope.org/metadata/fits?obs_id=%d -O %d.metafits",obsid,obsid);
         printf("Running : %s",szCmd);
         int ret_val = system(szCmd);      
      }
   }
   m_Metafits = new CBgFits( szMetaFitsPath );   
   if( !m_Metafits ){
      printf("ERROR : could not allocate object new for m_Metafits\n");
      exit(-1);
   }
   
   
   // CBgFits metafits( szMetaFitsPath );
   if( m_Metafits->ReadFits( szMetaFitsPath , 1 , 0 ) ){
       printf("ERROR : could not read metafits file %s - this error can not be tolerated -> exiting now !!!\n",szMetaFitsPath);       
       exit(-1);
   }

   // READ CHANNELS using metafits file :   
   HeaderRecord* pKeyword = m_Metafits->GetKeyword( "CHANNELS" );
   if( pKeyword ){
       CMyStrTable items;
       MyParser szPars = pKeyword->Value.c_str();
       szPars.GetItems(items,",&");
       printf("Keyword CHANNELS = %s -> split into %d items = |",pKeyword->Value.c_str(),((int)items.size()));
       if( items.size() > 0 ){
           m_CoarseChannels.clear();
           m_CoarseChannels.assign(items.size(),0);
           for(int i=0;i<items.size();i++){
               printf("%s|",items[i].c_str());
               m_CoarseChannels[i] = atol( items[i].c_str() );
           }       
           printf("\n");
           sort( m_CoarseChannels.begin() , m_CoarseChannels.end() );
                      
           m_FreqUpperMHz = MWA_COARSE_CHANNEL/2.00 + m_CoarseChannels[m_CoarseChannels.size()-1]*MWA_COARSE_CHANNEL;
           printf("m_CoarseChannels.size() =  %d -> m_FreqUpperMHz = %.2f MHz\n",((int)m_CoarseChannels.size()),m_FreqUpperMHz);
       }
   }else{
       printf("WARNING : could not get keyword CHANNELS in metafits file %s\n",szMetaFitsPath);
   }
   printf("INFO : inttime  = %.2f seconds\n",m_Metafits->inttime);
//   printf("INFO : unixtime = %d\n",m_Metafits->dtime_fs); 

   return true;
}

int CMWADataCube::FindReferenceSources( int border , const char* szOutDir )
{
   string mean_regfile = "mean.reg";
   string median_regfile = "median.reg";
//   const char* change_ext(const char* name,const char* new_ext,string& out);
   

   if( m_pMeanImage ){
      change_ext( m_pMeanImage->GetFileName(), "reg", mean_regfile );
      m_pMeanImage->FindSources( 5, border, mean_regfile.c_str() );           
   }else{
      printf("WARNING : m_pMeanImage = NULL -> cannot find reference sources in this image\n");
   }

   if( m_pMedianImage ){
      change_ext( m_pMedianImage->GetFileName(), "reg", median_regfile );
      m_pMedianImage->FindSources( 5, border, median_regfile.c_str() );
   }else{
      printf("WARNING : m_pMedianImage = NULL -> cannot find reference sources in this image\n");
   }
   
   return 1;
}


int CMWADataCube::get_freq_list( std::vector<double>& freq_list )
{
   freq_list.clear();
   for(int i=0;i<m_CoarseChannels.size();i++){
     double freq_mhz = m_CoarseChannels[i] * MWA_COARSE_CHANNEL;      
     
     freq_list.push_back( freq_mhz );
   }
   
   return freq_list.size();
}

// averages n_avg_timesteps timesteps around every single start time ;
// WARNING : overwrites the real value !!!
int CMWAFits::avg_in_time( int n_avg_timesteps , const char* szOutAvgFitsFile )
{
  if( n_avg_timesteps > 1 ){
     CMWAFits dyna_spec_tmp( "tmp", GetXSize(), GetYSize() );
     dyna_spec_tmp.SetValue( 0.00/0.00 );

     
     // loop over channels : 
     for(int y=0;y<GetYSize();y++){
        // loop over start time to average over specified number of steps (n_avg_timesteps) starting at given start time :
        for(int start_t=0;start_t<GetXSize();start_t++){
           double sum = 0;
           int count = 0;
           
           double start_value = getXY( start_t, y );           
           if( !isnan(start_value) ){
              for(int t=start_t;t<(start_t+n_avg_timesteps);t++){ 
                 if( t>=0 && t<GetXSize() ){
                    double value = getXY( t, y );

                    if( !isnan(value) ){                 
                       sum += value;
                       count++;
                    }
                 }
              }
           
//              if( count > 0 ){
//                 sum = sum / count;
//              }
              sum = sum / count; // allow NaNs in the gaps between the observations 
              dyna_spec_tmp.setXY( start_t, y, sum );
           }else{
              printf("DEBUG : NaN value pixel (%d,%d) ignored\n",start_t,y);
           }
        }
     }


     for(int y=0;y<GetYSize();y++){
        for(int x=0;x<GetXSize();x++){
           double new_value = dyna_spec_tmp.getXY( x, y );

/*           if( gZeros2NaNs ){
              if( new_value == 0.00 || fabs(new_value) < gZeroDistanceThreshold ){
                 new_value = (0.00/0.00); // works better than FP_NAN
              }           
           }*/

           setXY( x, y, new_value );
        }
     }
     printf("INFO : overwritten the original dynamic spectrum with one averaged over %d timesteps\n",n_avg_timesteps);


     if( szOutAvgFitsFile ){
        dyna_spec_tmp.SetKeys( GetKeys() );   
        dyna_spec_tmp.WriteFits( szOutAvgFitsFile );
        printf("INFO(upgrade 202009) : dynamic spectrum averaged over %d timesteps was writted to %s FITS file\n",n_avg_timesteps,szOutAvgFitsFile );
     }
   }  
   
   return n_avg_timesteps;  
}




void CMWAFits::CalcStatPerChannel()
{
   if( m_StdDevVsChannel.size() > 0 ){
      printf("WARNING : statistics per frequency channel has already been calculated -> exiting CMWAFits::CalcStatPerChannel now\n");
      return;
   }

   printf("DEBUG : CMWAFits::CalcStatPerChannel sizes = %d / %d\n",int(m_StdDevVsChannel.size()),int(m_MeanVsChannel.size()));
   m_StdDevVsChannel.assign( GetYSize(), 0.00/0.00 ); // initialise with NaNs    
   m_MeanVsChannel.assign( GetYSize(), 0.00/0.00 ); // initialise with NaNs    

   FILE* out_f = fopen("rmsmean_vs_channel.txt","w");
   fprintf(out_f,"# freq_channel MEAN   RMS\n");
   
   double sum_rms = 0.00, sum2_rms = 0.00;
   int nChannels = GetYSize();
   int count_rms = 0;
   
   for(int y=0;y<nChannels;y++){
      double sum = 0.00, sum2 = 0.00;
      
      int count = GetXSize();
      for(int x=0;x<count;x++){
         double val = getXY(x,y);
         
         sum += val;
         sum2 += (val*val);
      }            
      
      double mean = sum/count;
      double rms = sqrt( sum2/count - mean*mean );
      if( !isnan(rms) ){
         sum_rms += rms;
         sum2_rms += (rms*rms);     
         count_rms += 1;
      }
      
      m_StdDevVsChannel[y] = rms;
      m_MeanVsChannel[y] = mean;
      fprintf(out_f,"%d %.8f %.8f\n",y,mean,rms);
   }   
   
   fclose(out_f);
      
   m_MeanRMS = sum_rms / count_rms;
   m_RmsRMS = sqrt( sum2_rms/count_rms - m_MeanRMS*m_MeanRMS);
   
   printf("INFO : RMS for freq. channel = %.8f +/- %.8f\n",m_MeanRMS,m_RmsRMS);
}

int  CMWAFits::AutoFlagChannels( double threshold_sigma )
{
   if( m_FlaggedChannels.size() > 0 ){
      printf("WARNING : auto-flagging already done before\n");
      return 0;
   }

   CalcStatPerChannel();
   
   FILE* out_f = fopen("auto_flagged_channels.txt","w");
   double threshold = threshold_sigma*m_RmsRMS + m_MeanRMS;
   int nChannels = GetYSize();   
   for(int y=0;y<nChannels;y++){
      printf("DEBUG : ch = %d , stddev = %.8f , rms_snr = %.8f vs. threshold = %.2f\n",y,m_StdDevVsChannel[y],(m_StdDevVsChannel[y] - m_MeanRMS)/m_RmsRMS,threshold_sigma);
      if( isnan(m_StdDevVsChannel[y]) || fabs((m_StdDevVsChannel[y] - m_MeanRMS)/m_RmsRMS) > threshold_sigma ){
         printf("\tINFO : channel %d flagged\n",y);
         m_FlaggedChannels.push_back(y);
         fprintf(out_f,"%d\n",y);
      }      
   }   


   // 
   // minimum number of ok channels to reset the RFI block / band : 
   int min_ok_channels = 5;
   int start_rfi_ch = -1;
   int end_rfi_ch = -1;
   int ok_block=0;
   int bad_block=0;
   
   for(int y=0;y<nChannels;y++){
      bool is_flagged = ( find_value( m_FlaggedChannels, y ) >= 0 );
      
      if( is_flagged ){
         if( start_rfi_ch < 0 ){
            start_rfi_ch = y;
         }
         end_rfi_ch = y;
         bad_block++;
         
         if( ok_block > 0 ){
            if( ok_block < min_ok_channels ){
               printf("DEBUG : channel = %d , ok_block = %d -> flagging range %d - %d\n",y,ok_block,start_rfi_ch,y);
               // flag channels up to and including the current channel :
               for(int ch=start_rfi_ch;ch<=y;ch++){
                  if( find_value( m_FlaggedChannels, ch ) < 0 ){
                     m_FlaggedChannels.push_back( ch );
                     fprintf(out_f,"%d\n",ch);
                  }
               }
            }
            ok_block = 0; // reset OK block as the current channel is BAD
         }
      }else{
         ok_block++;
         if( ok_block >= min_ok_channels ){
            // ok block large enough to keep it unflagged
            start_rfi_ch = -1;
            end_rfi_ch   = -1;            
         }
      }
   }

   fclose(out_f);

   // fill bool table
   m_IsFlagged.assign( nChannels, false );
   for(int y=0;y<nChannels;y++){
      if( find_value(m_FlaggedChannels,y) >= 0 ){
         m_IsFlagged[y] = true;
      }
   }   
   
   return m_FlaggedChannels.size();
}

int CMWAFits::InitFlagsByRanges( std::vector<CFlagRanges>& bad_channel_ranges )
{
   m_IsFlagged.assign( GetYSize(), false );
   
   for(int r=0;r<bad_channel_ranges.size();r++){
      CFlagRanges& range = bad_channel_ranges[r];
      for(int ch=range.m_Start;ch<=range.m_End;ch++){
         m_IsFlagged[ch] = true;
      }
   }
   
   // setting zero values to flagged channels
   for(int y=0;y<GetYSize();y++){
      if( m_IsFlagged[y] ){
         for(int x=0;x<GetXSize();x++){
            setXY(x,y,0.00);
         }
      }
   }
   
   WriteFits("tmp.fits");
   
   return 1;
}

// is channel flagged :
bool CMWAFits::is_channel_flagged( int channel )
{
//   int ch_idx = find_value( m_FlaggedChannels, channel );
   if( channel >= 0 && channel < m_IsFlagged.size() ){   
      return m_IsFlagged[ channel ];
   }
   
   printf("ERROR : error in code, channel %d outside of range in the array %d - %d\n",channel,0,int(m_IsFlagged.size()));
   return false;
}


void CMWAFits::ReadCRValues()
{
   crval1 = 0.00;   // start physical value
   cdelt1 = 1.00; // delta of physical value
   crpix1 = 1.00; // which pixel is the start 
   crval2 = 0.00;   // start physical value
   cdelt2 = 1.00; // delta of physical value
   crpix2 = 1.00; // which pixel is the start 
   // HeaderRecord* CBgFits::GetKeyword(const char *keyword)
   HeaderRecord* pKey = NULL;
   pKey = GetKeyword( "CRVAL1" );
   if( pKey ){
     crval1 = atof( pKey->Value.c_str() );
   }
   pKey = GetKeyword( "CDELT1" );
   if( pKey ){
     cdelt1 = atof( pKey->Value.c_str() );
   }

   pKey = GetKeyword( "CRPIX1" );
   if( pKey ){
     crpix1 = atof( pKey->Value.c_str() );
   }

   pKey = GetKeyword( "CRVAL2" );
   if( pKey ){
     crval2 = atof( pKey->Value.c_str() );
   }

   pKey = GetKeyword( "CDELT2" );
   if( pKey ){
     cdelt2 = atof( pKey->Value.c_str() );
   }

   pKey = GetKeyword( "CRPIX2" );
   if( pKey ){
     crpix2 = atof( pKey->Value.c_str() );
   }
}

