#include <iostream>
using namespace std;

#include <stripes/DenseMatrix.h>
#include "Viewer.h"
using namespace DDG;

#include <unistd.h>

int main( int argc, char** argv )
{
   srand( time( NULL ));
   Viewer viewer;

#ifdef SP_COMPILE_COMMAND_LINE
   if( argc != 2 )
   {
      cerr << "usage: " << argv[0] << " input.obj" << endl;
      return 1;
   }

   viewer.mesh.read( argv[1] );
#endif

   viewer.init();

   return 0;
}

