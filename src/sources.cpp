#include "sources.h"
#include <math.h>
#include <stdlib.h>

void CSources::Add( int x, int y, double flux, double rms, double mean )
{
    CSource new_source( x, y, flux, rms, mean );
    push_back( new_source );
}

CSource* CSources::find( int x, int y, double radius )
{
   for(int i=0;i<size();i++){
      CSource& source = (*this)[i];
      
      double dist = sqrt( (source.x-x)*(source.x-x) + (source.y-y)*(source.y-y) );
      if (dist <= radius ){
         return &source;
      }
   }
   
   return NULL;
}

