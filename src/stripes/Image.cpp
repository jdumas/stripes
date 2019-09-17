#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

#include "Image.h"

namespace DDG
{
   Image :: Image( size_t width, size_t height )
   : w( width ), h( height ), pixels( w*h*3 )
   {}

   float& Image :: operator()( size_t x, size_t y )
   // accesses pixel (x,y)
   {
      return pixels[ x + y*w ];
   }

   const float& Image :: operator()( size_t x, size_t y ) const
   // accesses pixel (x,y)
   {
      return pixels[ x + y*w ];
   }

   float Image :: sample( float x, float y ) const
   // samples image at (x,y) using bilinear filtering
   {
      const Image& I( *this );
      float ax = x - floor( x );
      float ay = y - floor( y );
      float bx = 1.f - ax;
      float by = 1.f - ay;
      size_t x0 = (size_t) floor( x );
      size_t y0 = (size_t) floor( y );
      size_t x1 = x0 + 1;
      size_t y1 = y0 + 1;

      clamp( x0, y0 );
      clamp( x1, y1 );

      return by * ( bx * I(x0,y0) + ax * I(x1,y0) ) +
             ay * ( bx * I(x0,y1) + ax * I(x1,y1) ) ;
   }

   size_t Image :: width( void ) const
   // returns image width
   {
      return w;
   }

   size_t Image :: height( void ) const
   // returns image height
   {
      return h;
   }

   class TGAHeader
   // header format for Truevision TGA images
   {
      public:
         char  idFieldSize;
         char  colorMapType;
         char  dataTypeCode;
         short colorMapOrigin;
         short colorMapLength;
         char  colorMapEntrySize;
         short xOrigin;
         short yOrigin;
         short width;
         short height;
         char  bitsPerPixel;
         char  imageSpecification;
   };

   void Image :: read( const char* filename )
   // loads an image file in Truevision TGA format
   // (must be uncompressed RGB image with 24 or 32 bits per pixel)
   {
      ifstream in( filename, ios_base::binary );

      if( !in.is_open() )
      {
         cerr << "Error: could not open file " << filename << " for input!" << endl;
         exit( 1 );
      }

      // read header
      TGAHeader header;
      in.read( (char*) &(header.idFieldSize),        1 );
      in.read( (char*) &(header.colorMapType),       1 );
      in.read( (char*) &(header.dataTypeCode),       1 );
      in.read( (char*) &(header.colorMapOrigin),     2 );
      in.read( (char*) &(header.colorMapLength),     2 );
      in.read( (char*) &(header.colorMapEntrySize),  1 );
      in.read( (char*) &(header.xOrigin),            2 );
      in.read( (char*) &(header.yOrigin),            2 );
      in.read( (char*) &(header.width),              2 );
      in.read( (char*) &(header.height),             2 );
      in.read( (char*) &(header.bitsPerPixel),       1 );
      in.read( (char*) &(header.imageSpecification), 1 );

      w = header.width;
      h = header.height;

      // validate data type
      const char uncompressedRGB = 2;
      if( header.dataTypeCode != uncompressedRGB ||
          ( header.bitsPerPixel != 24 &&
            header.bitsPerPixel != 32 ))
      {
         cerr << "Error: input must be uncompressed RGB image with 24 or 32 bits per pixel." << endl;
         exit( 1 );
      }

      // read identification field (unused)
      vector<char> idField( header.idFieldSize );
      in.read( &idField[0], header.idFieldSize );

      // read color map data (unused)
      if( header.colorMapType == 1 )
      {
         int bytesPerEntry = header.colorMapEntrySize / 8;
         int colorMapSize = header.colorMapLength * bytesPerEntry;
         vector<char> colorMapData( colorMapSize );
         in.read( &colorMapData[0], colorMapSize );
      }

      // read pixel data
      size_t n = w*h*header.bitsPerPixel/8;
      vector<unsigned char> pixelData( n );
      in.read( (char*) &pixelData[0], n );

      // convert pixel data to floating point
      pixels.resize( n );
      for( int i = 0; i < n; i++ )
      {
         pixels[i] = (unsigned char) ((double) pixelData[i] / 255.);
      }

      if( header.bitsPerPixel == 24 )
      {
         for( size_t i = 0; i < n/3; i++ )
         {
            swap( pixels[i*3+0], pixels[i*3+2] );
         }
      }
   }

   void Image :: write( const char* filename ) const
   // writes an image file in Truevision TGA format
   // (uncompressed RGB image with 24 bits per pixel)
   {
      ofstream out( filename, ios_base::binary );

      if( !out.is_open() )
      {
         cerr << "Error: could not open file " << filename << " for output!" << endl;
         exit( 1 );
      }

      TGAHeader header;
      header.idFieldSize = 0;
      header.colorMapType = 0;
      header.dataTypeCode = 2;
      header.colorMapOrigin = 0;
      header.colorMapLength = 0;
      header.colorMapEntrySize = 0;
      header.xOrigin = 0;
      header.yOrigin = 0;
      header.width = (short) w;
      header.height = (short) h;
      header.bitsPerPixel = 24;
      header.imageSpecification = 0;

      // write header
      out.write( (char*) &(header.idFieldSize),        1 );
      out.write( (char*) &(header.colorMapType),       1 );
      out.write( (char*) &(header.dataTypeCode),       1 );
      out.write( (char*) &(header.colorMapOrigin),     2 );
      out.write( (char*) &(header.colorMapLength),     2 );
      out.write( (char*) &(header.colorMapEntrySize),  1 );
      out.write( (char*) &(header.xOrigin),            2 );
      out.write( (char*) &(header.yOrigin),            2 );
      out.write( (char*) &(header.width),              2 );
      out.write( (char*) &(header.height),             2 );
      out.write( (char*) &(header.bitsPerPixel),       1 );
      out.write( (char*) &(header.imageSpecification), 1 );

      // convert pixel data from floating point
      vector<unsigned char> pixelData( w*h*3 );
      for( int i = 0; i < w*h*3; i++ )
      {
         pixelData[i] = (unsigned char)( pixels[i] * 255. );
      }

      // write pixel data
      out.write( (char*) &pixelData[0], w*h*3 );
   }

   void Image :: clamp( size_t& x, size_t& y ) const
   // clamps coordinates to range [0,w-1] x [0,h-1]
   {
      x = std::max( 0ull, std::min( w-1, x ));
      y = std::max( 0ull, std::min( h-1, y ));
   }
}

