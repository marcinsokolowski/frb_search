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

enum eNormalisationType { eNoNorm=0, eNormMinMax=1, eNormRMS=2 };
eNormalisationType gNormaliseProfile = eNoNorm;

string gOutputDir="./";
string gInputFitsFile = "dyna_spec.fits";
string gOutFitsFile="series.fits";
double gPeakAvoidInterval = 0.02; // seconds - where to calculate MEDIAN and RMSIQR 
double gOffPulseInterval=0.02; // in seconds
long int gOffPulseIntervalInSamples=(long int)(0.02/((1.08*64*14)/1000.00)); // in ms 
int gLocalMedianSamples=40;
double gSamplingTime = (1.08*64*14)/1000000.00; // in seconds

bool gSaveFITS=false;

// testing DM
double gTestDM = 56.73;
double gStartDM = -1;
double gEndDM   = -1;
double gStepDM  = 0.001;

// metadata :
double gUNIXTIME=-1;

// debug level 
int gVerb=0;

void usage()
{
   printf("dedisperse_fits FITS_FILE_DYNA-SPEC\n");
   printf("-o OUTDIR : output directory [default %s]\n",gOutputDir.c_str());
   printf("-S START_DM [default %.6f]\n",gStartDM);
   printf("-E START_DM [default %.6f]\n",gEndDM);
   printf("-s STEP_DM [default %.6f]\n",gStepDM);
   printf("-D DM [default %.6f]\n",gTestDM);
   printf("-U UNIXTIME [default %.8f]\n",gUNIXTIME);

   exit(-1);
}


void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "a:ho:S:E:s:D:U:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      printf("DEBUG : opt = %c, optarg = %s\n",opt,optarg);
   
      switch (opt) {
         case 'o' :
            gOutputDir = optarg;
            printf("DEBUG : set output directory to %s\n",gOutputDir.c_str());
            break;
            
         case 'D' :
            gTestDM = atof( optarg );
            break;
            
         case 'S' :
            gStartDM = atof( optarg );
            break;
            
         case 'E' :
            gEndDM = atof( optarg );
            break;
            
         case 's' :
            gStepDM = atof( optarg );
            break;
            
         case 'U' :
            gUNIXTIME = atof( optarg );
            break;
            
         case 'h':
            usage();
            break;
         default:  
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }   
}
 

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("Input fits file  = %s\n",gInputFitsFile.c_str());
   printf("Output directory = %s\n",gOutputDir.c_str());  
   printf("DM range = %.6f - %.6f\n",gStartDM,gEndDM);
   printf("DM step = %.6f\n",gStepDM);
   printf("Test DM = %.6f\n",gTestDM);
   printf("#####################################\n");   
}

void dedisperse_fits( CBgFits& dynaspec, double DM )
{
   double delta_x = dynaspec.cdelt1*1000.00; // Multiplied by 100 to convert to ms
   double starting_frequency = 199.60214170; // Starting frequency used by PRESTO 
   double delta_y = dynaspec.cdelt2;

   std::vector<double> frequencies;
   int n_ch = dynaspec.GetYSize();
   for(int ch=0;ch<n_ch;ch++){
      frequencies.push_back(ch*delta_y);
   }


   double freq = starting_frequency + (n_ch-1)*delta_y;
   
   // de-dispersion :
   for(int ch=0;ch<n_ch;ch++){
      double fch1 = starting_frequency + (ch * delta_y);
      
      // tdelay = Decimal(Decimal(4.15)* 10**6 * DM * ((fch1**-2) - (freq**-2)))
      double tdelay = 4.15*1000000.00*DM*(1.00/(fch1*fch1) - 1.00/(freq*freq));
      
      double shift = (tdelay / delta_x);
      int shift_int = -int(round(shift));
      
      if( gVerb >= 1 ){
         printf("DEBUG : channel = %d -> shift = %.4f -> %d\n",ch,shift,shift_int);
      }
      dynaspec.roll( ch, shift_int );      
   }
}

#define REPLACE_ELEMS( tab, pos1, pos2 ) { tmp=tab[pos1]; tab[pos1]=tab[pos2]; tab[pos2]=tmp; } 

void my_sort_float( float* ftab, long cnt )
{
        float divider = ftab[0];

        int beg = 1;
        int end = cnt-1;
        float tmp;

        if(cnt){        
                while(beg<=end){
                        if(ftab[beg]>divider){
                                if(ftab[end]<=divider){
                                        REPLACE_ELEMS( ftab, beg, end )
                                        beg++;
                                        end--;
                                }else{
                                        end--;
                                }
                        }else{          
                                beg++;
                                if(ftab[end]>divider)
                                        end--;
                        }
                }
                if(end!=0){
                        REPLACE_ELEMS( ftab, end, 0)
                }

                my_sort_float( ftab, end );
                my_sort_float( &(ftab[beg]), cnt-beg );
        }

}


double get_trim_median( double n_sigma_iqr, double* tab, int& cnt, double& sigma_iqr )
{   
   double* newtab = new double[cnt];

   int q75= (int)(cnt*0.75);
   int q25= (int)(cnt*0.25);
   double iqr = tab[q75]-tab[q25];
   sigma_iqr = iqr/1.35;
   double range = sigma_iqr*n_sigma_iqr;
   double median = tab[(int)cnt/2];

   int newcnt=0;
   for(int i=0;i<cnt;i++){
      if( fabs(tab[i]-median) <= range ){
         newtab[newcnt] = tab[i];
         newcnt++;
      }
   }

   double ret=newtab[newcnt/2];

   // returning smaller array :
   cnt = newcnt;
   for(int i=0;i<newcnt;i++){
      tab[i] = newtab[i];
   }   

   delete [] newtab;
   return ret;
}

// INPUT  : 
// intab  : sorted table                                                 
// cnt    : number of elements in a table
// n_iter : number of iterations 
void GetAvgEstimator( float* intab, int cnt, int n_iter, double& sigma_iqr, double& median ) // , int x )
{
   double* tab = new double[cnt]; // working copy 
   for(int i=0;i<cnt;i++){
      tab[i] = intab[i];
   }

   int newcnt=cnt;
   median = tab[cnt/2];
   
   for(int i=0;i<n_iter;i++){
      median = get_trim_median( 3.00, tab, newcnt, sigma_iqr );
   }   

   delete [] tab;
}





void calc_median_rmsiqr( float* data_buffer, long int n_floats, long int start, long int end, double& median, double& rmsiqr )
{
   long int tmp_buffer_size = (end-start+1);
   float* tmp_buffer = new float[tmp_buffer_size];
   int tmp_count=0;
   
   for(long int i=start;i<end;i++){
      tmp_buffer[tmp_count] = data_buffer[i];
      tmp_count++;
   }
   
   my_sort_float( tmp_buffer, tmp_count );
   GetAvgEstimator( tmp_buffer, tmp_count, 5, rmsiqr, median );
   
   
   delete [] tmp_buffer;
}

void calc_median_rmsiqr( std::vector<double>& data_buffer, long int start, long int end, double& median, double& rmsiqr )
{
   int n_floats = data_buffer.size();
   long int tmp_buffer_size = (end-start+1);
   float* tmp_buffer = new float[tmp_buffer_size];
   int tmp_count=0;
   
   for(long int i=start;i<end;i++){
      tmp_buffer[tmp_count] = data_buffer[i];
      tmp_count++;
   }
   
   my_sort_float( tmp_buffer, tmp_count );
   GetAvgEstimator( tmp_buffer, tmp_count, 5, rmsiqr, median );
   
   
   delete [] tmp_buffer;
}

int find_peak( std::vector<double>& timeseries, int& min_t )
{
   double maxval=-1;
   int max_t=-1;
   
   min_t=-1;
   double minval = 1e20;
      
   for(int t=0;t<timeseries.size();t++){
      if( timeseries[t] > maxval ){
         maxval = timeseries[t];
         max_t = t;
      }
      if( timeseries[t] < minval ){
         minval = timeseries[t];
         min_t = t;
      }
   }
   
   return max_t;
}

double calc_snr_for_dm( CBgFits& infits, double dm, double external_rmsiqr=-1.00 )
{
  printf("\n\n");
  printf("Testing dm = %.6f\n",dm);

  std::vector<double> timeseries;
  timeseries.assign( infits.GetXSize(), 0 );

  dedisperse_fits( infits, dm );  
  
  for(int t=0;t<infits.GetXSize();t++){
     double sum = 0.00;
     for(int ch=0;ch<infits.GetYSize();ch++){
        sum += infits.getXY(t,ch);
     }
     
     timeseries[t] = sum;
  }
  
  int min_t=-1;
  int max_t = find_peak( timeseries, min_t );
  if ( max_t >= 0 ){
     printf("Peak found at time index = %d (value = %.8f)\n",max_t,timeseries[max_t]);
  }
  double minval = timeseries[min_t];
  double maxval = timeseries[max_t];
  
  if( gNormaliseProfile == eNormMinMax ){     
     for(int t=0;t<infits.GetXSize();t++){
        timeseries[t] = (timeseries[t]-minval)/(maxval-minval);
     }
  }
  
  // calculate MEDIAN/RMSIQR 
  double post_pulse_median=-1000, post_pulse_rmsiqr=-1000;
  long int post_pulse_start = max_t + gOffPulseIntervalInSamples;
  long int post_pulse_end = post_pulse_start + gLocalMedianSamples;
  calc_median_rmsiqr( timeseries, post_pulse_start, post_pulse_end, post_pulse_median, post_pulse_rmsiqr );
  
  double pre_pulse_median=-1000, pre_pulse_rmsiqr=-1000;
  long int pre_pulse_start = max_t - gOffPulseIntervalInSamples - gLocalMedianSamples;
  long int pre_pulse_end = pre_pulse_start + gLocalMedianSamples;
  calc_median_rmsiqr( timeseries, pre_pulse_start, pre_pulse_end, pre_pulse_median, pre_pulse_rmsiqr );
  
  double mean_median = (post_pulse_median + pre_pulse_median)/2.00;
  double mean_rmsiqr = (post_pulse_rmsiqr + pre_pulse_rmsiqr)/2.00;
  double diff_median = (post_pulse_median - pre_pulse_median);
  double diff_rmsiqr = (post_pulse_rmsiqr - pre_pulse_rmsiqr);
  double diff_median_frac = (post_pulse_median - pre_pulse_median)/mean_median;
  double diff_rmsiqr_frac = (post_pulse_rmsiqr - pre_pulse_rmsiqr)/mean_rmsiqr;
  printf("DEBUG : median/rmsiqr calculated %ld samples before and after the peak (%.6f sec)\n",gOffPulseIntervalInSamples,gSamplingTime*gOffPulseIntervalInSamples);
  printf("DEBUG : pre median/rmsiqr = %.8f / %.8f , post median/rmsiqr = %.8f / %.8f -> DIFF %.8f / %.8f ( %.8f / %.8f )\n",pre_pulse_median, pre_pulse_rmsiqr, post_pulse_median, post_pulse_rmsiqr, diff_median, diff_rmsiqr, diff_median_frac, diff_rmsiqr_frac );
  
  char szPulseOutFile[128],szPrePulseOutFile[128],szPostPulseOutFile[128],szSNROutFile[128],szMedianFile[128];
  sprintf(szSNROutFile,"%s/snr_vs_time_dm%.6f.txt",gOutputDir.c_str(),dm);
  sprintf(szPulseOutFile,"%s/timeseries_dm%.6f.txt",gOutputDir.c_str(),dm);
  sprintf(szPrePulseOutFile,"%s/prepulse_median_rms_dm%.6f.txt",gOutputDir.c_str(),dm);
  sprintf(szPostPulseOutFile,"%s/postpulse_median_rms_dm%.6f.txt",gOutputDir.c_str(),dm);
  sprintf(szMedianFile,"%s/median_rmsiqr.txt",gOutputDir.c_str());
  FILE* pulse_out_f = fopen(szPulseOutFile,"w");
  FILE* snr_out_f = fopen(szSNROutFile,"w");
  FILE* prepulse_out_f = fopen(szPrePulseOutFile,"w");
  FILE* postpulse_out_f = fopen(szPostPulseOutFile,"w");
  FILE* median_out_f = fopen(szMedianFile,"w");

  double rmsiqr_to_use = mean_rmsiqr;
  if( external_rmsiqr > 0 ){
     rmsiqr_to_use = external_rmsiqr;
  }
  
  double max_snr = -1;
  for(int t=0;t<timeseries.size();t++){
     double t_ms = t*gSamplingTime*1000.00;
     double t_sec = t*gSamplingTime;
     double snr = (timeseries[t]-mean_median)/rmsiqr_to_use;
     
     if( snr > max_snr ){
        max_snr = snr;
     }
  
     fprintf(pulse_out_f,"%.6f %.8f\n",t_sec,timeseries[t]);     
     fprintf(snr_out_f,"%.6f %.8f\n",t_sec,snr);
     fprintf(prepulse_out_f,"%.8f %.8f %.8f\n",t_sec,pre_pulse_median,pre_pulse_rmsiqr);
     fprintf(postpulse_out_f,"%.8f %.8f %.8f\n",t_sec,post_pulse_median,post_pulse_rmsiqr);
     fprintf(median_out_f,"%.6f %.8f %.8f\n",t_sec,mean_median,mean_rmsiqr);
  }
  fclose( median_out_f );
  fclose( pulse_out_f );
  fclose( snr_out_f );
  fclose( prepulse_out_f );
  fclose( postpulse_out_f );

/*  if( gNormaliseProfile == eNormMinMax && gNormaliseFits ){
     for(int y=0;y<infits.GetYSize();y++){
        for(int x=0;y<infits.GetXSize();x++){
           double val = infits.getXY(x,y);
           val = (val-minval)/(maxval-minval);
           infits.setXY(x,y,val);
        }
     }
  }*/

  if( gSaveFITS ){
     char gOutFullFitsPath[256];
     sprintf(gOutFullFitsPath,"%s/dm%.6f.fits",gOutputDir.c_str(),dm);  
     infits.WriteFits( gOutFullFitsPath );
  }

  

  return max_snr;
}

int main(int argc,char* argv[])
{
  time_t start_time = get_dttm();
  if( (argc<2 || ( argc>=2 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--h")==0))) ){
     usage();
  }

  gInputFitsFile = argv[1];  
  if( argc>=3 ){
     gOutFitsFile = argv[2];
  }

  parse_cmdline(argc,argv); // -1,+1
  print_parameters();
  
  // create output directory :
  MyFile::CreateDir( gOutputDir.c_str() );
  printf("INFO : created directory %s\n",gOutputDir.c_str());

  CBgFits infits;
  if( infits.ReadFits( gInputFitsFile.c_str() ) ){
     printf("ERROR : could not read FITS file %s\n",gInputFitsFile.c_str());
     exit(-1);
  }
  infits.ReadCRValues();  
  printf("OK : read input FITS file %s\n",gInputFitsFile.c_str());
  
  std::vector<double> timeseries;
  timeseries.assign( infits.GetXSize(), 0 );
  for(int t=0;t<infits.GetXSize();t++){
     double sum = 0.00;
     for(int ch=0;ch<infits.GetYSize();ch++){
        sum += infits.getXY(t,ch);
     }

     timeseries[t] = sum;
  }
  double median_start, rmsiqr_start;
  calc_median_rmsiqr( timeseries, 10, 10+gLocalMedianSamples, median_start, rmsiqr_start );
  printf("MEDIAN_START = %.8f, RMSIQR_START = %.8f\n",median_start,rmsiqr_start);


  char szOutFile[128];
  sprintf(szOutFile,"%s/maxsnr_vs_dm.txt",gOutputDir.c_str()); 
  FILE* outf = fopen(szOutFile,"w");

  if( gStartDM > 0.00 && gEndDM > 0.00 && gStepDM > 0.00 ){
     CBgFits tmpfits( infits.GetXSize(), infits.GetYSize() );
     if( tmpfits.ReadFits( gInputFitsFile.c_str() ) ){
        printf("ERROR : could not read FITS file %s\n",gInputFitsFile.c_str());
        exit(-1);
     }
     tmpfits.ReadCRValues();  

     int size = infits.GetXSize()*infits.GetYSize();
     float* tmp_data = tmpfits.get_data();

     double maxmaxsnr = -1, best_dm = -1;  
     double dm = gStartDM;
     while( dm <= gEndDM ){
        // WARNING : de-dispersion works in-place (changes FITS file);
        // copy from original image 
        memcpy( tmp_data, infits.get_data(), sizeof(float)*size );
        
        double maxsnr = calc_snr_for_dm( tmpfits, dm, rmsiqr_start );  
        fprintf(outf,"%.6f %.6f\n",dm,maxsnr);
        
        if( maxsnr > maxmaxsnr ){
           maxmaxsnr = maxsnr;
           best_dm = dm;
        }

        dm += gStepDM;
     }
     
     FILE* outf2 = fopen("DM.txt","w");
     fprintf(outf2,"%.8f %.8f\n",gUNIXTIME,best_dm);
     fclose(outf2);
  }else{
     if( gTestDM > 0 ){
        double maxsnr = calc_snr_for_dm( infits, gTestDM );  
        fprintf(outf,"%.6f %.6f\n",gTestDM,maxsnr);
     }
  }
  
  fclose(outf);
}
