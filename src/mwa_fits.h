
#ifndef _MWA_FITS_H__
#define _MWA_FITS_H__

#include <bg_fits.h>
#include "sources.h"

#include <string>
#include <vector>

extern int gDebugLocal;

enum eDumpType { eDumpType_Ptr=0, eDumpType_CenterPixel=1 };

class CPixel;
class CDedispSearch;

#define ZERO_VALUE -10000

class CFlagRanges
{
public:
   int m_Start;
   int m_End;
   
   CFlagRanges( int start=0, int end=0 )
    : m_Start(start), m_End(end)
   {
   }
};

class CMWAFits : public CBgFits
{
public :
   static bool m_bUseLocalThreshold;
    
   // static variables :
   // static double MWA_COARSE_CHANNEL; // = 1.28; # MHz 
   // static int    MWA_COARSE_CHANNELS;// = 24;   # Number of MWA coarse channels
   int m_ReadStatus;

   double m_Mean;
   double m_RMS;
   double m_Min;
   double m_Max;
   int m_UsedChannels;
   int m_Obsid;
   std::vector<string> m_OriginalImages;
   
   // flagged channels :
   std::vector<int> m_FlaggedChannels; // list of flagged channels
   std::vector<bool> m_IsFlagged; // table of all channels with true for flagged (bad) or false for ok channel
   // std::vector<CFlagRanges> m_FlagRanges;
   std::vector<double> m_StdDevVsChannel;
   std::vector<double> m_MeanVsChannel;
   double m_MeanRMS;
   double m_RmsRMS;
   
   
   // list of sources
   CSources m_Sources;
   
   static int m_Border;
   
   // start and end time:
   double m_StartTime;
   
   // WCS coordinates / CRVALUES 
   double crval1;   // CRVAL1
   double cdelt1;  // CDELT1
   double crpix1;   // CRPIX1
   double crval2;   // CRVAL2
   double cdelt2;  // CDELT2
   double crpix2;   // CRPIX2
   
   double get_x_physical_value( int x ){ double ret =  crval1 + double(x)*cdelt1; printf("%.8f + %.1f*%.8f = %.8f\n",crval1,double(x),cdelt1,ret); return ret; }
   double get_y_physical_value( int y ){ return crval2 + double(y)*cdelt2; }
   
   void ReadCRValues();

   CMWAFits( const char* fits_file=NULL, int xSize=0, int ySize=0 )
    : CBgFits( fits_file, xSize, ySize ), m_Mean(0.00), m_RMS(1e20), m_Min(0), m_Max(0), m_ReadStatus(-1), m_UsedChannels(0), m_Obsid(-1), m_MeanRMS(0.00), m_RmsRMS(1e20), m_StartTime(0.0),
      crval1(0.00), cdelt1(1.00), crpix1(1.00), crval2(0.00), cdelt2(1.00), crpix2(1.00)
   {
   }
   
   double GetStat( double& mean, double& rms, double& minval, double& maxval, CPixel& max_pixel, CPixel& min_pixel, 
                   CBgArray& max_along_x, vector<int>& max_along_x_x,
                   CBgArray& max_along_y, vector<int>& max_along_y_y,
                   int x_start=0, int y_start=0, int x_end=-1, int y_end=-1, 
                   double min_allowed_value=-1e20, double max_allowed_value=+1e20 );
   double GetStat( int border = 20 );
   int FindSources( double threshold_in_sigma, int border = 20, const char* szOutRegFile="sources.reg", bool bFluxInSigmas=true, int minX=0, int maxX=-1, 
                    CBgFits* pCountMap=NULL, CBgArray* p_rms_per_count=NULL, 
                    int max_count_param=-1 // pass how many MAX number pixels can be in dispersive delay sweep
                  );
   int FindSourcesSNR( double threshold_snr=5.00, int border = 20, const char* szOutRegFile="sources_snr.reg", int minX=0, int maxX=-1, CBgFits* pCountMap=NULL , bool bSavePhysical=false);
   
   
   // sets values smaller than threshold in sigmas to -1000 
   // returns number of pixels above threshold_in_sigma 
   int CutValues( double threshold_in_sigma /* was =2.00 */, vector<int>& good_columns, double zero_value = ZERO_VALUE, int rms_border=2 );
   
   // removes ZERO columns 
   int CopyWithoutBadColumns( CMWAFits& out_fits, vector<int>& good_columns );
   
   // average in time direction :
   int avg_in_time( int n_avg_timesteps , const char* szOutAvgFitsFile=NULL );

   // RFI per channel flagging    
   void CalcStatPerChannel();
   int AutoFlagChannels( double threshold_sigma=0.5 /* sigma */ );
   int InitFlagsByRanges( std::vector<CFlagRanges>& bad_channel_ranges );

   // is channel flagged :
   bool is_channel_flagged( int channel );
};


class CMWADataCube : public std::vector< std::vector< CMWAFits* > > 
{
public :
   // metafits :
   CBgFits* m_Metafits;  
   
   // mean image :
   CMWAFits* m_pMeanImage;
   
   // sum2 image 
   CMWAFits* m_pRMSImage;
   
   // median image :
   CMWAFits* m_pMedianImage;
   
   // flags for whole images :
   std::vector< std::vector< int > > m_ImageFlags;
   
   int m_Obsid;
   int m_Channels;
   int m_Timesteps;
   int m_StartTimeIndex;
   int m_FirstCorrectTimestep;
   int m_FirstCorrectChannel;
   int m_FirstCoarseChannel;
   vector<double> m_CoarseChannels; // TODO : should be double to be more general 
   vector<double> m_fCoarseChannelFreqMHz; // frequencies of channels as double , TODO : replace m_CoarseChannels with m_fCoarseChannelFreqMHz
   double m_FreqLowerMHz;
   double m_FreqUpperMHz;
   double m_FreqCenterMHz;
   double m_StartUnixTime;
   

   CMWADataCube();
   CMWADataCube( int obsid, int n_channels, int timesteps, int start_timeindex=0 );
   void Init( int obsid, int n_channels, int timesteps, int start_timeindex, double freq_center_mhz=-1, double freq_lower_mhz=-1, double freq_upper_mhz=-1 );
   ~CMWADataCube();
   
   CMWAFits* GetImage( int channel, int timeindex );
   int GetXSize(); // { return GetImage(0,0)->GetXSize(); }
   int GetYSize(); // { return GetImage(0,0)->GetYSize(); }
   void SetImage( int channel, int timeindex, CMWAFits* pBgFits ){
      ((*this)[channel][timeindex]) = pBgFits;
   }
   int GetDynamicSpectrum( int x, int y, CBgFits& out_dynaspec ); 

   void DebugDump( eDumpType dump_type = eDumpType_Ptr );   
   int Read( const char* dir_template="wsclean_timeindex%03d", 
             const char* image_template="wsclean_%d_timeindex%03d-%04d-I-dirty.fits",
             bool bExitOnReadError = true
           );
           
   // reading in channel (00000/ etc) and Time directory structures, it requires fits_list in each channel directory to exist;        
   int ReadChanTime( const char* dir_template="%05d",  // channel directory name
                     const char* image_template="dirty_image_%dT%d_real.fits",
                     bool bExitOnReadError = true
                   );
           
   int ReadFitsCube( const char* fits_list, double freq_lower_mhz, double freq_upper_mhz, 
                     int bAutoDetect=0, int bReadImage=1, int bIgnoreHeaderErrors=0 );
                     
   // save FITS files for invividual channels :                        
   int SaveChannelFitsFiles( int n_images=-1 );
           
   bool ReadMetaData( int obsid = -1 );
   bool GenerateMetaData( CBgFits& dyna_spec, CDedispSearch& dedisp, int obsid = -1, double start_freq=-1, double delta_freq=-1 );
             
   // automatically check number of available timesteps :
   static int count_time_steps( const char* dir_template="wsclean_timeindex%03d", int max_timesteps=1000  );
   static int DoesDirExist( const char* szDirPath );
             
   int FindReferenceSources( int border=20, const char* szOutDir="./" );
   
   int GetClosestImage( double unixtime, double& min_dist );
   int GetOptimalImage( double unixtime, double& min_dist );
   
   int get_freq_list( std::vector<double>& freq_list );
   
   // generating file files :
   int generate_filfile( int x=-1, int y=-1, const char* out_filename=NULL, double f_start=142.0760, double foff=-0.0040, double tstart=58323.2563, double bw=1.28, 
                         double inttime=0.010, bool bSaveRescaledFits=false );
   
   // Arithmetic operations :
   bool Operation( CMWADataCube& right , eCalcFitsAction_T oper );
   bool Avg( CMWADataCube& right );
};

extern double MWA_COARSE_CHANNEL;  //  = 1.28; // # MHz 
extern int    MWA_COARSE_CHANNELS; // = 24;  // # Number of MWA coarse channels 
extern int    gDebugLocal;

#endif