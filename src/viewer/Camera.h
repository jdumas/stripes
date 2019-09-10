// -----------------------------------------------------------------------------
// libDDG -- Camera.h
// -----------------------------------------------------------------------------
//
// Camera is used by Viewer to keep track of the view state; it also
// handles mouse input related to camera manipulation.
//

#ifndef DDG_CAMERA_H
#define DDG_CAMERA_H

#include <stripes/Quaternion.h>

#include <glad/glad.h>
#include <GL/glut.h>

namespace DDG
{
   class Camera
   {
      public:
         Camera( void );
         // constructor

         Quaternion clickToSphere( int x, int y );
         // projects a mous click onto the unit sphere

         void setView( void ) const;
         // applies the camera transformation to the OpenGL modelview stack

         void mouse( int button, int state, int x, int y );
         // handles mouse clicks

         void motion( int x, int y );
         // handles mouse drags

         void idle( void );
         // handles camera momentum

         void zoomIn( void );
         // moves viewer toward object

         void zoomOut( void );
         // moves viewer away from object

         void moveUp( void );
         void moveDown( void );
         void moveLeft( void );
         void moveRight( void );
         // translate in image plane

         Quaternion currentRotation( void ) const;
         // returns the rotation corresponding to the current mouse state

         Quaternion pClick;
         // mouse coordinates of current click

         Quaternion pDrag;
         // mouse coordinates of current drag

         Quaternion pLast;
         // mouse coordinates of previous drag

         Quaternion rLast;
         // previous camera rotation

         Quaternion momentum;
         // camera momentum

         int tLast;
         // time of previous drag

         double zoom, vZoom;
         // zoom and zoom velocity

         Vector eye;
         // camera position
   };
}

#endif

