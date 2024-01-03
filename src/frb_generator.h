#ifndef _FRB_GENERATOR_H__
#define _FRB_GENERATOR_H__

class CBgFits;

class CFrbGenerator 
{
public :
    static void AddSource( CBgFits& fits, int x, int y, double mean_flux, double int_start_time, double inttime );
    static void AddSourceTest( CBgFits& fits, int x, int y );
    static Double_t FRB_Pulse( Double_t* x, Double_t* y );
    static double GetMeanFluxDensity(  double freq_mhz, double start_time, double end_time, double int_start_time, double inttime, int bScattering , Double_t scattering_tau );
    static double calc_channel_start_time(  int coarse_channel, 
                          double freq1_mhz = 200.32, /* MHz - end of upper band of frequency range = 156*1.28 + 1.28/2.00 */
                          double dm=114.00           // DM for FRB171020 
                        );
                        
    static double FluxIntegral( double int_start_time , double int_end_time, double freq_mhz, double start_time,  double_t scattering_tau, double timestep = 1e-18 );

};    


#endif

