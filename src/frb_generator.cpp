// taken from pulse_injector/main.cpp
// g++ frb_generator.cpp -c -I ../fitslib/ -I$SRCDIR/cmn/baselib/ -D_UNIX

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

#ifdef _HAS_ROOT_
// root section
#include <TFile.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TNtuple.h> 
#include <TRandom.h>
#include <TGraph2DErrors.h>
#include <TVirtualFitter.h>
#include <TStopwatch.h>
#include <TError.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TGraph.h>
#else

#define Double_t double
#include <math.h>

#endif

#include "frb_generator.h"

#define N_CHANNELS 24

double gTimeIndexStepSec=0.5; // time interval of timeindex (default 0.5 seconds images)
int    gInsertCoarseChannel=17;
double gFreqEnd = 200.32; // changed on 20180830 from 200.00 MHz -> 200.32 MHz so that it is arrival time at the upper end of the coarse channel band 
                          // upper end of highest coarse channel ( = 156 *1.28 + 0.64 MHz )
                          // was 200.32 but now the same to agree with test : python ./dispersion2.py 114 200 170
double gDM      = 114.00;

double gTimeOffset=0.00;

int gDebug=1;

double gSpectralIndex = -1.8; // spectral index 

double gScatteringTime=-1.00;
int    bScattering=0;

void CFrbGenerator::AddSource( CBgFits& fits, int x, int y, double mean_flux, double int_start_time, double inttime )
{
   double current_value = fits.getXY( x, y );
   double new_value = current_value + mean_flux;
   double fluence = mean_flux * inttime;
   fits.setXY( x, y, new_value );
   printf("\t%.3f sec : Adding fluence = %.3f Jy sec -> mean flux = %.3f Jy sec / %.2f sec = %.4f Jy -> changing %.4f Jy -> %.4f Jy\n",(double)int_start_time,fluence,fluence,inttime,mean_flux,current_value,new_value);
}

void CFrbGenerator::AddSourceTest( CBgFits& fits, int x, int y )
{
   fits.setXY( x, y, 10 );
   fits.setXY( x, y+1, 5 );
   fits.setXY( x, y-1, 5 );
   fits.setXY( x-1, y, 5 );
   fits.setXY( x+1, y, 5 );
}

// No scattering tau = 0.050
// Scattering tau = 2 seconds ?
Double_t CFrbGenerator::FRB_Pulse( Double_t* x, Double_t* y )
{
   Double_t t = x[0];
   Double_t freq_mhz = y[1];
   Double_t freq_ghz = freq_mhz / 1000.00;

   Double_t t0 = y[2];           // start of the pulse

   // FRB171020 , 200 Jy ms @ 1.4 GHz 
   //             1500 Jy ms @ 185 MHz and spectral index  = -1
   //             11400 Jy ms @ 185 MHz and spectral index = -2
   //             7600 Jy ms @ 185 MHz and spectral index  = -1.8
   Double_t fluence_askap = 200;  // Jy x sec assuming spectral index  = -1.8
   Double_t spectral_index = gSpectralIndex; // was -1.8 
   Double_t fluence_Jy_ms = fluence_askap*pow( (freq_ghz / 1.400) , spectral_index ); // Jy ms 
   Double_t fluence = fluence_Jy_ms / 1000.00;  // Jy x second 
   // calculate amplitude based on frequency and ASKAP value at 1.4 MHz 
   // see /home/msok/Desktop/MWA/logbook/FRB_ASKAP_RW/paper/20180525_FRB_fluence_at_MWA.odt
   // Spectral index = -1.80 : min_fluence=6659.75 [Jy] @ 199.68 MHz, max_fluence=8874.62 [Jy] @ 170.24 MHz, central = 7596.61 [Jy] @ 185.60 
   // python ./frb_fluence.py 200 170.24 1400.00 114 500  | grep Spec   
//   if( gDebug > 0 ){
//       printf("\t\tDEBUG : t0 = %.2f sec , fluence( %.4f MHz ) = %.4f Jy ms  = %.4f Jy second\n",t0,freq_mhz,fluence_Jy_ms,fluence);
//   }
   

   Double_t tau = y[0] * pow( (200.00/freq_mhz) , 4 );      // 50ms - WARNING, when <0.050 - the integral value is incorrect !!! scattering time as a function of frequnecy using index = 4 
//   printf("DEBUG : tau(%.2f MHz) = %.4f seconds\n",freq_mhz,tau);
   Double_t Amplitude = fluence / tau; // assuming pulse : Amplitude x exp( - t / tau )

   if( t < t0 ){
     return 0;
   }

   double pulse = 0.00;
   if( t > t0 ){
       pulse = Amplitude* exp( - (t-t0) / tau );
   }
   return pulse;
}

double CFrbGenerator::FluxIntegral( double int_start_time , double int_end_time, double freq_mhz, double start_time,  double_t scattering_tau, double timestep /* = 1e-18 */ )
{
   double params[3];
   params[0] = scattering_tau;
   params[1] = freq_mhz;
   params[2] = start_time;;
   
   
   double sum = 0.00;
   double t = int_start_time;
   while( t < int_end_time ){
      double flux = FRB_Pulse( &t , params );
      sum += flux;
   
      t+= timestep;
   }
   
   double fluence = sum*timestep;
   
   return fluence;
}

double CFrbGenerator::GetMeanFluxDensity(  double freq_mhz, double start_time, double end_time, double int_start_time, double inttime, int bScattering , Double_t scattering_tau )
{
   Double_t sigma_disp     = 0.130; // 130 mJy (see /home/msok/Desktop/MWA/logbook/FRB_ASKAP_RW/paper/proposal/submission/ApJL_submission/Referee_comments/Averages_tests/referee_comments_on_avarege_vs_scattering
   Double_t sigma_half_sec = 0.580; // 580 mJy
   Double_t sigma_10dedisp = 0.092; //  92 mJy
   Double_t sigma_60dedisp = 0.084; //  84 mJy


#ifdef _HAS_ROOT_
   TF1* line = NULL; 
   if( bScattering <= 0 ){
       line = new TF1("FRB_Pulse2",FRB_Pulse,0,end_time,3);
       // tau without scattering :
       line->SetParameter(0,0.050); // 50ms - WARNING, when <0.050 - the integral value is incorrect !!!
   }else{
       line = new TF1("FRB_Pulse2",FRB_Pulse,0,end_time + scattering_tau*5.00 , 3);
       // tau with scattering :
       line->SetParameter(0, scattering_tau); 
       
   }
   line->SetParameter(1,freq_mhz);
   line->SetParameter(2,start_time); 

   // double expected_sigma_long_integration = sigma_disp / sqrt( n_averaged );
   // double observed_sigma_long_integration = interpolate( n_averaged, 10, sigma_10dedisp, 60, sigma_60dedisp );
   // double sigma_long_integration = expected_sigma_long_integration;
   // if( bDrawOptimal <= 0 ){
   //   sigma_long_integration = observed_sigma_long_integration;
   // }
   // double sigma_half_second   = sigma_disp;

   double fluence_half_second = line->Integral( int_start_time, int_start_time + inttime, 1e-18 );
   delete line;
#else
   double fluence_half_second = FluxIntegral( int_start_time, int_start_time + inttime, freq_mhz, start_time, scattering_tau, 1e-18 );
#endif   
   
   double mean_flux_density   = fluence_half_second / inttime;
   
   
   
   return mean_flux_density;
}

double CFrbGenerator::calc_channel_start_time( 
                          int coarse_channel, 
                          double freq1_mhz, /*= 200.32  MHz - end of upper band of frequency range = 156*1.28 + 1.28/2.00 */
                          double dm /*=114.00 */          // DM for FRB171020 
                        )
{
    int n_coarse = 23 - coarse_channel;
    double freq2_mhz = freq1_mhz - n_coarse*1.28;

    double freq2_ghz = freq2_mhz / 1000.00;
    double freq1_ghz = freq1_mhz / 1000.00;
    
    double delta_t_ms = 4.15 * dm * ( 1.00/(freq2_ghz*freq2_ghz) - 1.00/(freq1_ghz*freq1_ghz) );
    double delta_t_sec = delta_t_ms / 1000.00;
    
        
    return delta_t_sec;
}

/*
int main(int argc,char* argv[])
{
  // print_cmdline(argc,argv);
  if( argc>=2 && strncmp(argv[1],"-h",2)==0 ){
     usage();
  }  

  string dir_list_file = argv[1];
  // parse command line :
  parse_cmdline(argc-1,argv+1);
  print_parameters();
  
  vector<string> dir_list;
  int n_dirs = bg_read_list( dir_list_file.c_str(), dir_list );
  printf("Read directory list file %s with %d directories\n",dir_list_file.c_str(),n_dirs);

  

  // images :  
  vector< vector<CBgFits*> > data_list;
  data_list.reserve( n_dirs );
  vector<CBgFits*> empty_vector;
  
  
  for( int i=0; i<dir_list.size(); i++ ){
     const char* list_file = dir_list[i].c_str();
     // add new empty vector to fill with list of 24 coarse channel images for a given time index :
     data_list.push_back( empty_vector );
     
     
     mystring szFitsBase,szDrv,dirname,szExt;
     mystring szTmpPath=list_file;
     szTmpPath += ".txt"; // just artifical for next line with splitpath to work OK !!!
     szTmpPath.splitpath(szDrv,dirname,szFitsBase,szExt);         
     
     // data_list[i].reserve( N_CHANNELS );
     printf("Reading list file : %s (dirname = %s)\n",list_file,dirname.c_str());
     vector<string> fits_list;
     int n_images = bg_read_list( list_file , fits_list );
     printf("\tReading %d images form list %s -> dirname = %s\n",n_images,list_file,dirname.c_str());
     if( n_images != N_CHANNELS ){
         printf("ERROR : wrong number of images %d != expected %d (coarse channels)\n",n_images,N_CHANNELS);
         exit(-1);
     }
     
     for (int image=0;image<fits_list.size();image++){
         string full_path = dirname.c_str();
         full_path = full_path + "/" + fits_list[image];
     
         CBgFits* pFits = new CBgFits( full_path.c_str() );
         if( pFits-> ReadFits( full_path.c_str(), 1 ) ){
             printf("ERROR : error reading fits file %s\n",full_path.c_str() );
             exit(-1);
         }
         data_list[i].push_back( pFits );
         printf("\tFile %s added to list\n",full_path.c_str() );
     }                   
  }
  
  double gFrbtime = 0.00; // frb starts at the first image :
  int    source_x  = -1;
  int    source_y  = -1;
  
  // adding FRB signal for all coarse channels going down from highest frequency channel down to lowest frequency coarse channel :
  for(int cc=(N_CHANNELS-1);cc>=0;cc--){
      int n_coarse = 23 - cc;
      double cc_start_freq = gFreqEnd - n_coarse*1.28;
      double freq_mhz = 133*1.28 + cc*1.28; // 133 is start of cc=145 band 
      double time_shift = gTimeOffset + calc_channel_start_time( cc,  gFreqEnd, gDM );
      
      printf("Channel %d (upper end at %.2f MHz) : time = %.2f mili-seconds\n",cc,cc_start_freq,time_shift*1000.00);
      // double GetMeanFluxDensity(  double start_time, double end_time, double int_start_time, double inttime, int bScattering , Double_t scattering_tau )
      for(int timeindex=0;timeindex<data_list.size();timeindex++){
         double int_start_time = timeindex * gTimeIndexStepSec;
         double frb_mean_flux = GetMeanFluxDensity( freq_mhz, time_shift, time_shift + 50, int_start_time, gTimeIndexStepSec, bScattering, gScatteringTime );
         
         vector<CBgFits*>& coarse_channel_images = data_list[timeindex];
         CBgFits* pFitsToAddPulse = coarse_channel_images[cc];
         
         if( source_x < 0 || source_y < 0 ){
             source_x = pFitsToAddPulse->GetXSize()/2.00;
             source_y = pFitsToAddPulse->GetYSize()/2.00;
             printf("Automatically injected pulse position = (%d,%d)\n",source_x,source_y);             
         }
         AddSource( *pFitsToAddPulse , source_x, source_y, frb_mean_flux, int_start_time, gTimeIndexStepSec );
         pFitsToAddPulse->UpdateImage( pFitsToAddPulse->GetFileName(), pFitsToAddPulse->GetFileName() );
      }
  }
}
*/