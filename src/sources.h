#ifndef _SOURCES_H__
#define _SOURCES_H__

#include <vector>

class CSource
{
public :
   CSource( int _x, int _y, double _flux, double _rms=0.00, double _mean=0.00 ) 
   : x(_x), y(_y), flux(_flux), rms(_rms), mean(_mean)
   {}
   
   CSource( const CSource& right )
   {
      (*this) = right;         
   }
   
   CSource& operator=( const CSource& right )
   {
       x = right.x;
       y = right.y;
       flux = right.flux;
       rms = right.rms;
       mean = right.mean;
       
       return (*this);
   }
   
   void Update( int _x, int _y, double _flux ){
      x = _x;
      y = _y;
      flux = _flux;
   }

   double x;
   double y;
   double flux;
   
   double ra_deg;
   double dec_deg;
   
   double dm;
   double rms;
   double mean;
};

class CSources : public std::vector<CSource>
{
public :
   void Add( int x, int y, double flux, double rms=0.00, double mean=0.00 );
   CSource* find( int x, int y, double radius=3 );
};

#endif
