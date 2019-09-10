#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

#include "Viewer.h"
#include "Image.h"
#include "Complex.h"
#include "Utility.h"
#include "Image.h"
using namespace DDGConstants;

namespace DDG
{
   // declare static member variables
   Mesh Viewer::mesh;
   Viewer::RenderMode Viewer::mode = renderShaded;
   GLuint Viewer::surfaceDL = 0;
   int Viewer::windowSize[2] = { 512, 512 };
   Camera Viewer::camera;
   Shader Viewer::stripeShader, Viewer::smoothShader;
   bool Viewer::showSingularities = false;
   bool Viewer::showDirectionField = false;
   bool Viewer::animate = false;
   GLuint Viewer::texture = 0;
   bool Viewer::doPick = false;
   int Viewer::pickX = 0;
   int Viewer::pickY = 0;
   int Viewer::pickIndex = -1;
   bool Viewer::doneSolved = false;
   bool Viewer::isTrivialSection = false;

   void Viewer :: init( void )
   {
      restoreViewerState();
      initGLUT();
      initGL();
      initGLSL();
#ifndef SP_COMPILE_COMMAND_LINE
      mesh.read( DATA_FOLDER "bunny.obj" );
#endif
      updateDisplayList();
      glutMainLoop();
   }

   void Viewer :: initGL( void )
   {
      glClearColor( 1., 1., 1., 1. );
   }

   void Viewer :: initGLUT( void )
   {
      int argc = 0;
      vector< vector<char> > argv(1);

      if(gladLoadGL()) {
         // you need an OpenGL context before loading glad
         printf("I did load GL with no context!\n");
         exit(-1);
      }

      // initialize window
      glutInitWindowSize( windowSize[0], windowSize[1] );
      glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
      glutInit( &argc, (char**)&argv );
      glutCreateWindow( "StripePatterns" );

      // specify callbacks
      glutDisplayFunc  ( Viewer::display  );
      glutIdleFunc     ( Viewer::idle     );
      glutKeyboardFunc ( Viewer::keyboard );
      glutSpecialFunc  ( Viewer::special  );
      glutMouseFunc    ( Viewer::mouse    );
      glutMotionFunc   ( Viewer::motion   );

      // initialize menus
      int viewMenu = glutCreateMenu( Viewer::view );
      glutSetMenu( viewMenu );
      glutAddMenuEntry( "[s] Smooth Shaded",  menuSmoothShaded );
      glutAddMenuEntry( "[f] Wireframe",      menuWireframe    );
      glutAddMenuEntry( "[↑] Zoom In",        menuZoomIn       );
      glutAddMenuEntry( "[↓] Zoom Out",       menuZoomOut      );

      int mainMenu = glutCreateMenu( Viewer::menu );
      glutSetMenu( mainMenu );
      glutAddMenuEntry( "[space] Process Mesh", menuProcess    );
      glutAddMenuEntry( "[e] Edit Field",       menuTrivial    );
      glutAddMenuEntry( "[r] Reset Mesh",       menuResetMesh  );
      glutAddMenuEntry( "[w] Write Mesh",       menuWriteMesh  );
      glutAddMenuEntry( "[\\] Screenshot",      menuScreenshot );
      glutAddMenuEntry( "[esc] Exit",           menuExit       );
      glutAddSubMenu( "View", viewMenu );
      glutAttachMenu( GLUT_RIGHT_BUTTON );

      if(!gladLoadGL()) {
         printf("Something went wrong!\n");
         exit(-1);
      }
   }

   void Viewer :: initGLSL( void )
   {
      smoothShader.loadVertex( SHADER_FOLDER "vertex.glsl" );
      smoothShader.loadFragment( SHADER_FOLDER "smooth.glsl" );

      stripeShader.loadVertex( SHADER_FOLDER "vertex.glsl" );
      stripeShader.loadFragment( SHADER_FOLDER "stripe.glsl" );
   }

   void Viewer :: menu( int value )
   {
      switch( value )
      {
         case( menuProcess ):
            mProcess();
            break;
         case( menuTrivial ):
            mConvertToTrivialSection();
            break;
         case( menuResetMesh ):
            mResetMesh();
            break;
         case( menuWriteMesh ):
            mWriteMesh();
            break;
         case( menuScreenshot ):
            mScreenshot();
            break;
         case( menuExit ):
            mExit();
            break;
         default:
            break;
      }
   }

   void Viewer :: view( int value )
   {
      switch( value )
      {
         case( menuSmoothShaded ):
            mSmoothShaded();
            break;
         case( menuWireframe ):
            mWireframe();
            break;
         case( menuZoomIn ):
            mZoomIn();
            break;
         case( menuZoomOut ):
            mZoomOut();
            break;
         default:
            break;
      }
   }

   void Viewer :: mProcess( void )
   {
      updateStripePattern();
      updateDisplayList();
   }

   void Viewer :: mResetMesh( void )
   {
      mesh.reload();
      doneSolved = false;
      updateDisplayList();
   }

   void Viewer :: mWriteMesh( void )
   {
      mesh.write( "output.obj" );
      //writeSingularities();
   }

   void Viewer :: mExit( void )
   {
      storeViewerState();
      exit( 0 );
   }

   void Viewer :: mSmoothShaded( void )
   {
      mode = renderShaded;
      updateDisplayList();
   }

   void Viewer :: mWireframe( void )
   {
      mode = renderWireframe;
      updateDisplayList();
   }

   void Viewer :: mZoomIn( void )
   {
      camera.zoomIn();
   }

   void Viewer :: mZoomOut( void )
   {
      camera.zoomOut();
   }

   void Viewer :: mScreenshot( void )
   {
      static int index = 0;

      // get window width and height
      GLint view[4];
      glGetIntegerv( GL_VIEWPORT, view );
      int w = view[2];
      int h = view[3];

      // get pixels
      Image image( w, h );
      glReadPixels( 0, 0, w, h, GL_BGR, GL_FLOAT, &image(0,0) );

      stringstream filename;
      filename << "frames/viewer" << setw(8) << setfill( '0' ) << index << ".tga";
      image.write( filename.str().c_str() );

      index++;
   }

   void Viewer :: keyboard( unsigned char c, int x, int y )
   {
      switch( c )
      {
         case 's':
            mSmoothShaded();
            break;
         case 'f':
            mWireframe();
            break;
         case 'd':
            showDirectionField = !showDirectionField;
            updateDisplayList();
            break;
         case 'w':
            mWriteMesh();
            break;
         case 'e':
            mConvertToTrivialSection();
            break;
         case 'r':
            mResetMesh();
            break;
         case '\\':
            mScreenshot();
            break;
         case ' ':
            mProcess();
            break;
         case 'c':
            mesh.computeCurvatureAlignedSection();
            mesh.parameterize();
            updateDisplayList();
            doneSolved = true;
            showSingularities = true;
            showDirectionField = false;
            isTrivialSection = false;
            break;
         case 9: // tab
            animate = !animate;
            break;
         case '*':
            showSingularities = !showSingularities;
            updateDisplayList();
            break;
         case '(':
            rotateFieldLeft();
            break;
         case ')':
            rotateFieldRight();
            break;
         case '-':
         case '_':
            decreaseStripeFrequency();
            break;
         case '+':
         case '=':
            increaseStripeFrequency();
            break;
         case 27:
            mExit();
            break;
         case '1':
            mesh.nCoordinateFunctions = 1;
            break;
         case '2':
            mesh.nCoordinateFunctions = 2;
            break;
         default:
            break;
      }
   }

   void Viewer :: rotateStripes( void )
   {
      double theta = .25*M_PI/360.;
      Complex r( cos(theta), sin(theta) );

      for( VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++ )
      {
         v->directionField *= r;
      }

      mesh.computeParameterization( 0 );

      updateDisplayList();
   }

   void Viewer :: increaseStripeFrequency( void )
   {
      mesh.lambda *= 1.1;
      mesh.parameterize();
      updateDisplayList();
   }

   void Viewer :: decreaseStripeFrequency( void )
   {
      mesh.lambda /= 1.1;
      mesh.parameterize();
      updateDisplayList();
   }

   void Viewer :: special( int i, int x, int y )
   {
      bool shift = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);

      switch( i )
      {
         case GLUT_KEY_UP:
            if( shift ) { camera.moveUp(); }
            else        { camera.zoomIn(); }
            break;
         case GLUT_KEY_DOWN:
            if( shift ) { camera.moveDown(); }
            else        { camera.zoomOut(); }
            break;
         case GLUT_KEY_LEFT:
            if( shift ) { camera.moveLeft(); }
            break;
         case GLUT_KEY_RIGHT:
            if( shift ) { camera.moveRight(); }
            break;
         case 27:
            mExit();
            break;
         default:
            break;
      }
   }

   void Viewer :: display( void )
   {
      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


      glMatrixMode( GL_PROJECTION );
      glLoadIdentity();

      GLint viewport[4];
      glGetIntegerv( GL_VIEWPORT, viewport );
      double aspect = (double) viewport[2] / (double) viewport[3];
      const double clipNear = .01;
      const double clipFar = 10000.;

      gluPerspective( 45., aspect, clipNear, clipFar );

      glMatrixMode( GL_MODELVIEW );
      glLoadIdentity();

      Quaternion    eye = Vector( 0., 0., -2.5*camera.zoom );
      Quaternion center = Vector( 0., 0., 0. );
      Quaternion     up = Vector( 0., 1., 0. );

      gluLookAt(    eye[1],    eye[2],    eye[3],
                 center[1], center[2], center[3],
                     up[1],     up[2],     up[3] );


      Quaternion r = camera.currentRotation();

      eye = r.bar() * eye * r;
      {
         smoothShader.enable();
         GLint uniformEye = glGetUniformLocation( smoothShader, "eye" );
         glUniform3f( uniformEye, eye[1], eye[2], eye[3] );
         smoothShader.disable();
      }
      {
         stripeShader.enable();
         GLint uniformEye = glGetUniformLocation( stripeShader, "eye" );
         glUniform3f( uniformEye, eye[1], eye[2], eye[3] );
         stripeShader.disable();
      }

      Quaternion light = Vector( -1., 1., -2. );
      light = r.bar() * light * r;
      {
         smoothShader.enable();
         GLint uniformLight = glGetUniformLocation( smoothShader, "light" );
         glUniform3f( uniformLight, light[1], light[2], light[3] );
         smoothShader.disable();
      }
      {
         stripeShader.enable();
         GLint uniformLight = glGetUniformLocation( stripeShader, "light" );
         glUniform3f( uniformLight, light[1], light[2], light[3] );
         stripeShader.disable();
      }

      camera.setView();

      stripeShader.enable();

      if( !doPick )
      {
         drawSurface();
         glutSwapBuffers();
      }

      stripeShader.disable();

      if( doPick )
      {
         pick();
      }

   }

   void Viewer :: drawSurface( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glEnable( GL_DEPTH_TEST );
      glEnable( GL_LIGHTING );

      glCallList( surfaceDL );
      drawSelection();

      glPopAttrib();
   }

   void Viewer :: drawMesh( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glEnable( GL_POLYGON_OFFSET_FILL );
      glPolygonOffset( 1., 1. );

      GLint uniformStripeAmount = glGetUniformLocation( stripeShader, "stripeAmount" );

      glUniform1f( uniformStripeAmount, showDirectionField?0.:1. );
      drawPolygons();
      glUniform1f( uniformStripeAmount, 1. );

      stripeShader.disable();
      smoothShader.enable();
      drawFieldSingularities();
      smoothShader.disable();

      glDisable( GL_POLYGON_OFFSET_FILL );

      if( mode == renderWireframe )
      {
         stripeShader.disable();
         drawWireframe();
      }

      drawDirectionField();

      glPopAttrib();
   }

   void Viewer :: drawPolygons( void )
   {
      glColor3f( 1., .4, 0. );
      GLint uniformStripeTex = glGetUniformLocation( stripeShader, "stripeTex" );
      glUseProgram( stripeShader );
      glUniform1i( uniformStripeTex, 0 );
      glActiveTexture( GL_TEXTURE0 );
      glBindTexture( GL_TEXTURE_2D, texture );
      glUseProgram( stripeShader );

      glBegin( GL_TRIANGLES );
      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         if( f->isBoundary() ) continue;

         double k = f->fieldIndex(2.);
         double nu = f->paramIndex[0];
         double nv = f->paramIndex[1];

         glMultiTexCoord1d( GL_TEXTURE2, k );

         if( k == 0. || !doneSolved )
         {
            if( mode == renderWireframe )
            {
               Vector N = f->normal();
               glNormal3dv( &N[0] );
            }

               int i = 0;
               HalfEdgeCIter he = f->he;
               do
               {
                  if( mode != renderWireframe )
                  {
                     Vector N = he->vertex->normal();
                     glNormal3dv( &N[0] );
                  }

                  Complex g = he->texcoord;
                  glTexCoord4d( g.re, g.im, nu, nv );
                  glMultiTexCoord3d( GL_TEXTURE1, i==0, i==1, i==2 );
                  glMultiTexCoord2dv( GL_TEXTURE3, &he->vertex->parameterization.re );

                  Vector p = he->vertex->position;
                  glVertex3dv( &p[0] );

                  i++;
                  he = he->next;
               }
               while( he != f->he );
         }
         else // singular triangle
         {
            // Get the three half edges.
            HalfEdgeCIter hij = f->he;
            HalfEdgeCIter hjk = hij->next;
            HalfEdgeCIter hkl = hjk->next;

            // Get the three vertices.
            VertexCIter vi = hij->vertex;
            VertexCIter vj = hjk->vertex;
            VertexCIter vk = hkl->vertex;

            // Get the three parameter values---for clarity, let "l"
            // denote the other point in the same fiber as "i".  Purely
            // for clarity, we will explicitly define the value of psi
            // at l, which of course is always just the conjugate of the
            // value of psi at i.
            Complex psiI = vi->parameterization;
            Complex psiJ = vj->parameterization;
            Complex psiK = vk->parameterization;
            Complex psiL = psiI.bar();

            double cIJ = ( hij->edge->he != hij ? -1. : 1. );
            double cJK = ( hjk->edge->he != hjk ? -1. : 1. );
            double cKL = ( hkl->edge->he != hkl ? -1. : 1. );

            // Get the three omegas, which were used to define our energy.
            double omegaIJ = hij->omega();
            double omegaJK = hjk->omega();
            double omegaKL = hkl->omega();

            // Here's the trickiest part.  If the canonical orientation of
            // this last edge is from l to k (rather than from k to l)...
            omegaKL *= cKL;
            // SIMPLER // if( cKL == -1. )
            // SIMPLER // {
            // SIMPLER //    // ...then the value of omega needs to be negated, since the
            // SIMPLER //    // original value we computed represents transport away from
            // SIMPLER //    // vertex i rather than the corresponding vertex l
            // SIMPLER //    omegaKL = -omegaKL;
            // SIMPLER // }
            // Otherwise we're ok, because the original value was computed
            // starting at k, which is exactly where we want to start anyway.

            // Now we just get consecutive values along the curve from i to j to k to l.
            // (The following logic was already described in our routine for finding
            // zeros of the parameterization.)
            if( hij->crossesSheets() )
            {
               psiJ = psiJ.bar();
               omegaIJ =  cIJ * omegaIJ;
               omegaJK = -cJK * omegaJK;
            }

            // Note that the flag hkl->crossesSheets() is the opposite of what we want here:
            // based on the way it was originally computed, it flags whether the vectors at
            // Xk and Xi have a negative dot product.  But here, we instead want to know if
            // the vectors at Xk and Xl have a negative dot product.  (And since Xi=-Xl, this
            // flag will be reversed.)
            if( !hkl->crossesSheets() )
            {
               psiK = psiK.bar();
               omegaKL = -cKL * omegaKL;
               omegaJK =  cJK * omegaJK;
            }

            // From here, everthing gets computed as usual.
            Complex rij( cos(omegaIJ), sin(omegaIJ) );
            Complex rjk( cos(omegaJK), sin(omegaJK) );
            Complex rkl( cos(omegaKL), sin(omegaKL) );

            double sigmaIJ = omegaIJ - ((rij*psiI)/psiJ).arg();
            double sigmaJK = omegaJK - ((rjk*psiJ)/psiK).arg();
            double sigmaKL = omegaKL - ((rkl*psiK)/psiL).arg();
            //double xi = sigmaIJ + sigmaJK + sigmaKL;

            double betaI = psiI.arg();
            double betaJ = betaI + sigmaIJ;
            double betaK = betaJ + sigmaJK;
            double betaL = betaK + sigmaKL;
            double betaM = betaI + (betaL-betaI)/2.;

            // // to verify, one checks that these values are both integers (they are!)
            // cout << betaM/M_PI << endl;
            // cout << (betaI+betaL)/M_PI << endl;

            Vector pi = vi->position;
            Vector pj = vj->position;
            Vector pk = vk->position;
            Vector pm = ( pi + pj + pk ) / 3.;

            Vector Ni = vi->normal();
            Vector Nj = vj->normal();
            Vector Nk = vk->normal();
            Vector Nm = f->normal();

            glNormal3dv( &Ni.x ); glMultiTexCoord3d( GL_TEXTURE1, 1., 0., 0. ); glTexCoord4d( betaI, 0., 0., 0. ); glVertex3dv( &pi.x );
            glNormal3dv( &Nj.x ); glMultiTexCoord3d( GL_TEXTURE1, 0., 1., 0. ); glTexCoord4d( betaJ, 0., 0., 0. ); glVertex3dv( &pj.x );
            glNormal3dv( &Nm.x ); glMultiTexCoord3d( GL_TEXTURE1, 0., 0., 0. ); glTexCoord4d( betaM, 0., 0., 0. ); glVertex3dv( &pm.x );

            glNormal3dv( &Nj.x ); glMultiTexCoord3d( GL_TEXTURE1, 0., 1., 0. ); glTexCoord4d( betaJ, 0., 0., 0. ); glVertex3dv( &pj.x );
            glNormal3dv( &Nk.x ); glMultiTexCoord3d( GL_TEXTURE1, 0., 0., 1. ); glTexCoord4d( betaK, 0., 0., 0. ); glVertex3dv( &pk.x );
            glNormal3dv( &Nm.x ); glMultiTexCoord3d( GL_TEXTURE1, 0., 0., 0. ); glTexCoord4d( betaM, 0., 0., 0. ); glVertex3dv( &pm.x );

            glNormal3dv( &Nk.x ); glMultiTexCoord3d( GL_TEXTURE1, 0., 0., 1. ); glTexCoord4d( betaK, 0., 0., 0. ); glVertex3dv( &pk.x );
            glNormal3dv( &Ni.x ); glMultiTexCoord3d( GL_TEXTURE1, 1., 0., 0. ); glTexCoord4d( betaL, 0., 0., 0. ); glVertex3dv( &pi.x );
            glNormal3dv( &Nm.x ); glMultiTexCoord3d( GL_TEXTURE1, 0., 0., 0. ); glTexCoord4d( betaM, 0., 0., 0. ); glVertex3dv( &pm.x );
         }
      }
      glEnd();

      glDisable( GL_TEXTURE_2D );
   }

   void Viewer :: drawSelection( void )
   {
      int nF = mesh.faces.size();
      if( pickIndex < 0 || pickIndex >= nF ) return;
      FaceCIter f = mesh.faces[pickIndex].he->face;

      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glDisable( GL_LIGHTING );
      glEnable( GL_DEPTH_TEST );
      glColor4f( 0., 0., 1., .25 );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      glLineWidth( 8. );

      stripeShader.disable();

      glBegin( GL_TRIANGLES );
      HalfEdgeCIter he = f->he;
      do
      {
         glVertex3dv( &he->vertex->position.x );
         he = he->next;
      }
      while( he != f->he );
      glEnd();

      glPopAttrib();
   }


   void Viewer :: drawWireframe( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glDisable( GL_LIGHTING );
      glColor4f( 0., 0., 0., .5 );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      glDepthMask( 0 );

      glBegin( GL_LINES );
      for( EdgeCIter e  = mesh.edges.begin();
                     e != mesh.edges.end();
                     e ++ )
      {
         glVertex3dv( &e->he->vertex->position[0] );
         glVertex3dv( &e->he->flip->vertex->position[0] );
      }
      glEnd();

      glPopAttrib();
   }

   void Viewer :: drawDirectionField( void )
   {
      if( !showDirectionField ) return;

      glPushAttrib( GL_ALL_ATTRIB_BITS );

      stripeShader.disable();

      glDisable( GL_LIGHTING );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      glEnable( GL_LINE_SMOOTH );
      glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
      glColor4f( 0., 0., 1., .5 );
      glLineWidth( 2. );

      double meanEdgeLength = 0.;
      for( EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++ )
      {
         meanEdgeLength += e->length();
      }
      meanEdgeLength /= (double) mesh.edges.size();

      double k = mesh.fieldDegree;

      glBegin( GL_LINES );
      for( VertexCIter v = mesh.vertices.begin();
           v != mesh.vertices.end();
           v ++ )
      {
         Vector p = v->position;
         double r = .75 * meanEdgeLength;

         for( double n = 0; n < k; n++ )
         {
            Vector Z = v->fieldVector( k, n ).unit();
            Vector q = p + r*Z;
            glVertex3dv( &p[0] );
            glVertex3dv( &q[0] );
         }
      }
      glEnd();

      glPopAttrib();
   }

   void Viewer :: drawFieldSingularities( void )
   {
      if( !showSingularities ) return;

      glPushAttrib( GL_ALL_ATTRIB_BITS );
      glMatrixMode( GL_MODELVIEW );

      for( FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++ )
      {
         double n = f->fieldIndex( mesh.fieldDegree );
         if( n != 0. )
         {
            if( n < 0. ) glColor3f( 0., 0., 1. );
            if( n > 0. ) glColor3f( 1., 0., 0. );

            Vector b = f->barycenter();
            const double radius = .02;
            const double slices = 24;
            const double stacks = 48;
            glPushMatrix();
            glTranslatef( b.x, b.y, b.z );
            glutSolidSphere( radius, slices, stacks );
            glPopMatrix();
         }
      }

      glPopAttrib();
   }

   void Viewer :: drawIsolatedVertices( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      // draw with big, round, red dots
      glPointSize( 5 );
      glHint( GL_POINT_SMOOTH_HINT, GL_NICEST );
      glEnable( GL_POINT_SMOOTH );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      glColor4f( 1., 0., 0., 1. ); // red

      glBegin( GL_POINTS );
      for( VertexCIter v  = mesh.vertices.begin();
                       v != mesh.vertices.end();
                       v ++ )
      {
         if( v->isIsolated() )
         {
            glVertex3dv( &v->position[0] );
         }
      }
      glEnd();

      glPopAttrib();
   }

   void Viewer :: updateDisplayList( void )
   {
      if( surfaceDL )
      {
         glDeleteLists( surfaceDL, 1 );
         surfaceDL = 0;
      }

      surfaceDL = glGenLists( 1 );

      glNewList( surfaceDL, GL_COMPILE );
      drawMesh();
      glEndList();
   }

   void Viewer :: mouse( int button, int state, int x, int y )
   {
      camera.mouse( button, state, x, y );

      if( (glutGetModifiers()&GLUT_ACTIVE_ALT) && state == GLUT_UP )
      {
         doPick = true;
         pickX = x;
         pickY = y;
         display();
      }
   }

   void Viewer :: motion( int x, int y )
   {
      camera.motion( x, y );
   }

   void Viewer :: idle( void )
   {
      camera.idle();
      if( animate )
      {
         //scaleStripes();
         rotateStripes();
      }
      glutPostRedisplay();
   }

   void Viewer :: storeViewerState( void )
   {
      ofstream out( ".viewer_state.txt" );

      out << camera.rLast[0] << endl;
      out << camera.rLast[1] << endl;
      out << camera.rLast[2] << endl;
      out << camera.rLast[3] << endl;
      out << camera.eye.x << endl;
      out << camera.eye.y << endl;
      out << camera.zoom << endl;

      GLint view[4];
      glGetIntegerv( GL_VIEWPORT, view );
      out << view[2] << endl;
      out << view[3] << endl;

      out << (int) mode << endl;
   }

   void Viewer :: restoreViewerState( void )
   {
      ifstream in( ".viewer_state.txt" );

      if( !in.is_open() ) return;

      in >> camera.rLast[0];
      in >> camera.rLast[1];
      in >> camera.rLast[2];
      in >> camera.rLast[3];
      in >> camera.eye.x;
      in >> camera.eye.y;
      in >> camera.zoom;

      in >> windowSize[0];
      in >> windowSize[1];

      int m;
      in >> m;
      mode = (RenderMode) m;
   }

   void Viewer :: updateStripePattern( void )
   {
      // mesh.computeCurvatureAlignedSection();
      mesh.computeSmoothestSection();

      // mesh.extractSingularities();
      // mesh.computeTrivialSection();
      // mesh.alignTrivialSection();

      mesh.parameterize();

      cerr << "Drawing mesh..." << endl;

      doneSolved = true;
      showSingularities = true;
      showDirectionField = false;
      isTrivialSection = false;
   }

   void Viewer :: createTexture( void )
   {
      // load an image
      Image image;
      image.read( "data/GridRedBlue.tga" );
      int N = image.width();

      if( texture )
      {
         glDeleteTextures( 1, &texture );
         texture = 0;
      }

      glGenTextures( 1, &texture );
      glBindTexture( GL_TEXTURE_2D, texture );
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
      glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
      gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA, N, N, GL_RGB, GL_FLOAT, &image(0,0) );
   }

   void Viewer :: pick( void )
   {
      unsigned char c[3];
      unsigned char& R( c[0] );
      unsigned char& G( c[1] );
      unsigned char& B( c[2] );

      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
      glPushAttrib( GL_ALL_ATTRIB_BITS );
      glEnable( GL_DEPTH_TEST );
      glBegin( GL_TRIANGLES );
      for( int i = 0; i < mesh.faces.size(); i++ )
      {
         R = ( i & 0xFF );
         G = ( i & 0xFF00 ) >> 8;
         B = ( i & 0xFF0000 ) >> 16;
         glColor3ub( R, G, B );
         HalfEdgeCIter he = mesh.faces[i].he;
         do
         {
            glVertex3dv( &he->vertex->position.x );
            he = he->next;
         }
         while( he != mesh.faces[i].he );
      }
      glEnd();
      glPopAttrib();

      GLint viewport[4];
      glGetIntegerv( GL_VIEWPORT, viewport );
      GLint h = viewport[3];

      glReadPixels( pickX, h-pickY, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, c );
      int oldPickIndex = pickIndex;
      pickIndex = R | (G<<8) | (B<<16);

      static int nPicks = 0;
      nPicks++;
      if( nPicks%2 == 0 )
      {
         if( isTrivialSection )
         {
            updateSingularities( oldPickIndex, pickIndex );
         }
      }

      //printSelectionData();

      doPick = false;
      display();
   }

   void Viewer :: printSelectionData( void )
   {
      int nF = mesh.faces.size();
      if( pickIndex < 0 || pickIndex >= nF ) return;
      FaceCIter f = mesh.faces[pickIndex].he->face;
      HalfEdgeCIter hij = f->he;
      HalfEdgeCIter hjk = hij->next;
      HalfEdgeCIter hkl = hjk->next;

      cerr << "FACE " << pickIndex << endl;
      cerr << "===============================" << endl;
      cerr << "field index: " << f->fieldIndex(2.) << endl;
      cerr << "param index: " << f->paramIndex[0] << endl;
      cerr << "crossing: " << hij->edge->crossesSheets << " " << hjk->edge->crossesSheets << " " << hkl->edge->crossesSheets << endl;
      cerr << " flipped: " << (hij->edge->he!=hij) << " " << (hjk->edge->he!=hjk) << " " << (hkl->edge->he!=hkl) << endl;
      cerr << endl;
   }

   void Viewer :: writeSingularities( void )
   {
      ofstream outPos( "singularitiesPositive.obj" );
      ofstream outNeg( "singularitiesNegative.obj" );

      for( FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++ )
      {
         Vector c = f->barycenter();
         double k = f->fieldIndex(2.);
         if( k > 0. ) outPos << "v " << c.x << " " << c.y << " " << c.z << endl;
         if( k < 0. ) outNeg << "v " << c.x << " " << c.y << " " << c.z << endl;
      }
   }

   void Viewer :: mConvertToTrivialSection( void )
   {
      // Grab the singularities from the current field, and
      // construct the corresponding trivial connection---
      // from there the singularities can be edited.
      if( !isTrivialSection )
      {
         mesh.extractSingularities();
         mesh.computeTrivialSection();
         mesh.alignTrivialSection();
         mesh.extractSingularities();
         mesh.parameterize();
         isTrivialSection = true;
         updateDisplayList();
      }
   }

   void Viewer :: updateSingularities( int a, int b )
   // modify singularities for "trivial connection" direction field at the last two selected triangles
   {
      double& Ia( mesh.faces[a].singularIndex );
      double& Ib( mesh.faces[b].singularIndex );

      // if the singularities aren't distinct, do nothing
      if( a == b ) return;

      // if neither face is singular, make an equal and opposite pair
      if( Ia == 0. && Ib == 0. )
      {
         Ia =  .5;
         Ib = -.5;
      }

      // if the singularities have equal and opposite sign, cancel them
      else if( Ib == -Ia )
      {
         Ia = 0.;
         Ib = 0.;
      }

      // if one face is singular and the other is not, just move the singularity
      else if( (Ia == 0. && Ib != 0.) || (Ia != 0. && Ib == 0.) )
      {
         swap( Ia, Ib );
      }

      // otherwise, both faces are singular but have the same index---do nothing
      else
      {
         return;
      }

      // update the field and parameterization
      mesh.computeTrivialSection();
      mesh.extractSingularities();
      mesh.parameterize();
      updateDisplayList();
   }

   void Viewer :: rotateFieldLeft( void )
   {
      double theta = M_PI/32.;
      Complex z( cos(theta), sin(theta) );

      for( VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++ )
      {
         v->directionField *= z;
      }

      mesh.parameterize();
      updateDisplayList();
   }

   void Viewer :: rotateFieldRight( void )
   {
      double theta = -M_PI/32.;
      Complex z( cos(theta), sin(theta) );

      for( VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++ )
      {
         v->directionField *= z;
      }

      mesh.parameterize();
      updateDisplayList();
   }
}

      // ofstream out( "branchPoints.svg" );
      // out << "<svg>" << endl;
      // out << "<circle fill=\"none\" stroke=\"#000000\" cx=\"0.0\" cy=\"0.0\" r=\"0.5\"/>" << endl;
      // for( FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++ )
      // {
      //    if( f->normal().z < 0. ) continue;
      //    double nu = f->paramIndex[0];
      //    Vector p = f->barycenter();
      //    if( nu > 0 )
      //    {
      //       out << "<circle fill=\"#FF0000\" cx=\"" << p.x << "\" cy=\"" << p.y << "\" r=\"0.01\"/>" << endl;
      //    }
      //    if( nu < 0 )
      //    {
      //       out << "<circle fill=\"#0000FF\" cx=\"" << p.x << "\" cy=\"" << p.y << "\" r=\"0.01\"/>" << endl;
      //    }
      // }
      // out << "</svg>" << endl;
