#ifndef _DEDISP_SEARCH_H__
#define _DEDISP_SEARCH_H__

#include <map>
#include <vector>

#include "mwa_fits.h"
#include "pixel.h"

#define LOG_FILE_HEADER "# START_TIMESTEP   X   Y   Flux[Jy]    DM    SNR   PixelList\n"

class MyOFile;

/*class CPixel
{
public :
   int time;
   int cc;
   float flux;
   
   CPixel( int _time=0, int _cc=0, float _flux=0.00 ): 
   time(_time), cc(_cc), flux(_flux)
   {       
   }
};*/

class CDispersedPixels : public vector<CPixel> 
{
public :
   CDispersedPixels();
   bool add_unique( CPixel& new_pixel );
   void save_reg_file( const char* filename="dispersed_path.reg" );
   
   void save_text_file( const char* filename="dispersed_path.txt" );
   
//   CDispersedPixels& operator==( const CDispersedPixels& right ){
//      return ((vector<CPixel>&)(*this)) == ((vector<CPixel>&)right);
//   }
};

enum eAlgorithmType_T {eAlgoTypeBruteForce=0, eAlgoTypeTemplates=1 };

class CEvent
{
public :
    int start_timeidx;
    int x,y;
    double dm;
    double snr;
    double flux;
    eAlgorithmType_T algo_type;
    vector<CPixel> m_Pixels;
    
    CEvent( int _start_timeidx, int _x, int _y, eAlgorithmType_T _algo_type, double _dm=0.00, double _snr=0.00, double _flux=0.00 ) :
       start_timeidx( _start_timeidx ), x(_x), y(_y), algo_type( _algo_type ), dm(_dm), snr(_snr), flux(_flux)  
    {
    }
    
};

class CEventList : public std::vector<CEvent> 
{
public : 
    void add( int start_timeidx, int x, int y, eAlgorithmType_T _algo_type=eAlgoTypeBruteForce, double _dm=0.00, double _snr=0.00, double _flux=0.00 );   
};


class CDedispSearch 
{
public :
   CMWADataCube m_MWADataCube;
   
   // static member variables / parameters :
   static bool   m_bAllowNegativeDM;
   static double m_MinDM;
   static double m_MaxDM;
   static double m_StepDM;
   static double m_SNRThreshold;
   static double m_ThresholdToCutInSigmas;
   
   static double m_MinRMSOnSingle;
   static double m_MaxRMSOnSingle;
   static int    m_MinUsedChannels;
   static int    m_MinUsedImagesInDedisp;
   
   static int    m_SkipTimeSteps; // how many timesteps to skip 
   
   static bool   m_bAllowPreStartArrivals; // calculate DTS with negative arrival times (max negative value = max(dipersive delay))

   // saving results / logfiles 
   static string m_ResultsPath;
   
   static int m_bRunBruteForceAlgo;
   static int m_bRunDynSpecAlgo;
   
   // 
   static int m_bRejectRefSources;
   
   // If overwrite de-dispersed files 
   static bool m_bForceOverwrite;
   
   
   // debug dedispersion (log the files used to calculated de-dispersed image):
   static bool m_bDebugDedispersion;
   static double m_DebugDM;

   // static functions :
   static double dispersion_milisec( double freq1_mhz, double freq2_mhz, double dm );
   static double dispersion_freq_mhz( double freq1_mhz, double time_value_ms, double dm );
   static int find_value_double( vector<double>& tab, double value, double error=0.000001 );

   // temporary objects to store dedispersed images:
   vector<CMWAFits*> m_DedispersedImages;   
   int m_DedispersedImagesIndex;
   
   // templates of dedispersion steps in the dynamic spectrum :
   // keys are DMs in format XX.X (1 digit)
   std::map<string, CDispersedPixels > m_DispersionTemplates;
   double m_MaxDispersiveSweep;
   
   // event candidates :
   CEventList m_Events;
   
   CDedispSearch( int obsid, int n_channels, int timesteps, int start_timeindex=0 );
   ~CDedispSearch();


   void CleanDedispersed( int bDelete=0 );
   void CalcReferenceImages();
   void Dedisperse_OLD( double dm );  
   void Dedisperse( double dm );
   int FindTransientsDynSpec( double threshold_in_sigma, int border, int max_radius = 1, int min_good_pixels = 12 ); // use at least 12 channels 
   int AnalyseDynSpec( CMWAFits& dynspec, int x, int y, double threshold_in_sigma, MyOFile* candidates_file=NULL, int max_radius = 1, int min_good_pixels = 12, std::vector<double>* p_channels_list=NULL  );
   int FindTransients( double dm, double threshold_in_sigma=10, int border=20 );
   int Run( double threshold_in_sigma=10, int border=20, int min_good_pixels=12  );
   
   // 20191020 - get de-dispersed series, see LWA paper :
   int GetDedispersedSeries( CMWAFits& dynspec, CMWAFits& dedispersed_series, CMWAFits& out_count_map, CMWAFits& rms_map_new, std::vector<double>* p_channels_list=NULL, double threshold_in_sigma=5.00, 
                             MyOFile* candidates_file=NULL, int max_radius=1, double LogBelowDM=-1e6, bool bIgnoreMaxValue=false, bool bDoWeigthByRMS=true );   
   int ShowSweepPixels( CMWAFits& dynspec, int start_timeindex, MyOFile* candidates_file=NULL, std::vector<double>* p_channels_list=NULL );
   
   // templates calculations :
   int CalcTemplates( double inttime = -1.00, double freq_upper_end_mhz = -1 , std::vector<double>* p_channels_list=NULL );
   int CalcTemplates( CMWAFits& dynamic_spectrum, double inttime = -1, std::vector<double>* p_channels_list = NULL );
   int CalcTemplatesExact( CMWAFits& dynamic_spectrum );


   void PrintTemplates();
   int CalcTemplates( CMWADataCube& datacube, double inttime = -1.00 , double freq_upper_end_mhz = -1 , std::vector<double>* p_channels_list=NULL );
   
   // simulations :
   double InjectFRB( CMWAFits& dynspec, int start_timeindex=0, double value=10, bool overwrite=false, double snr=-100, int* end_timeindex=NULL );
};


#endif
