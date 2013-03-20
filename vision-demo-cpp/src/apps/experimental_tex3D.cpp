
#include <iostream>
#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <X11/Intrinsic.h>    /* Display, Window */
#include <GL/glx.h>           /* GLXContext */


using namespace std;

static GLuint  explodeProgram; // The geometry shader 'exploder' program
static GLint   locMVP;         // ModelViewProjection matrix uniform
static GLint   locMV;          // ModelView matrix uniform
static GLint   locNM;          // Normal matrix uniform

static GLint   szX=32;
static GLint   szY=32;
static GLint   szZ=32;

double _persp_m[16], _ortho_m[16];	// projection matrices

void  chernobylPointGrid( ) {


//  glEnable(GL_TEXTURE_3D);
//  glEnable(GL_FRAGMENT_PROGRAM_ARB);
//  glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, _prog[0]);
//  glActiveTextureARB(GL_TEXTURE0_ARB);
//  glBindTexture(GL_TEXTURE_3D, _txt[0]);

  glPointSize(5.0f);
  glColor3f(1.0f,0.5f,0.5f);
  glBegin(GL_POINTS); //starts drawing of points
    for( int i=0; i<szY; i++ ) {
      float y=i*2.0/(szY-1)-1.0f;
      for( int j=0; j<szX; j++ ) {
        float x=j*2.0/(szX-1)-1.0f;
        for( int k=0; k<szZ; k++ ) {
          float z=k*2.0/(szZ-1)-1.0f;
          glVertex3f(x,y,z);
          glTexCoord3f( (x+1.0)/2.0, (y+1.0)/2.0, (z+1.0)/2.0 );
        }
      }
    }
  glEnd();


  // if drawing something else, disable texture 3D first!
//  glDisable(GL_TEXTURE_3D);
//  glDisable(GL_FRAGMENT_PROGRAM_ARB);


}

unsigned int _txt[3];			  	// texture handles
unsigned int _prog[3];				// shader program handles
unsigned char*  _texture_data = NULL;
void Setup_RenderingContext(void)
{
  _texture_data = new unsigned char [szX*szY*szZ*4];
  int index =0;
  for( int i=0; i<szY; i++ ) {
    float y=i*2.0/(szY-1)-1.0f;
    for( int j=0; j<szX; j++ ) {
      float x=j*2.0/(szX-1)-1.0f;
      for( int k=0; k<szZ; k++ ) {
        float z=k*2.0/(szZ-1)-1.0f;
        _texture_data[index+0]=255*(x<0); // B
        _texture_data[index+1]=255*(y>0); // G
        _texture_data[index+2]=255*(z>0); // R
        _texture_data[index+3]=32; // A
        index+=4;
            //(j*spectrum_width*tex_steps+i*spectrum_width+k)*4;
      }
    }
  }

  // Background
  glClearColor(0.0f, 0.0f, 0.0f, 0.5f );
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_3D);
  glEnable(GL_BLEND);

  glGenTextures(2, &_txt[0]);

  glActiveTextureARB(GL_TEXTURE0_ARB);
  glBindTexture(GL_TEXTURE_3D, _txt[0]);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S,GL_CLAMP);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T,GL_CLAMP);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R,GL_CLAMP);

  glActiveTextureARB(GL_TEXTURE1_ARB);
  glBindTexture(GL_TEXTURE_3D, _txt[1]);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S,GL_CLAMP);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T,GL_REPEAT);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R,GL_REPEAT);

  glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, szX,
               szY, szZ, 0, GL_RGBA, GL_UNSIGNED_BYTE, _texture_data);

//  shaderManager.InitializeStockShaders();
//  viewFrame.MoveForward(4.0f);

//  // Make the torus
//  gltMakeTorus(torusBatch, .70f /*R1*/, 0.10f/*r2*/, 11, 7);

//  explodeProgram = gltLoadShaderTripletWithAttributes("GSExplode.vs",
//                                                      "GSExplode.gs",
//                                                      "GSExplode.fs",
//                                                      2,
//                                                      GLT_ATTRIBUTE_VERTEX, "vVertex",
//                                                      GLT_ATTRIBUTE_NORMAL, "vNormal");

//  locMVP = glGetUniformLocation(explodeProgram, "mvpMatrix");
//  locMV  = glGetUniformLocation(explodeProgram, "mvMatrix");
//  locNM  = glGetUniformLocation(explodeProgram, "normalMatrix");
//  locPushOut = glGetUniformLocation(explodeProgram, "push_out");
}

// Cleanup
void Shutdown_RenderingContext(void)
{
  glDeleteProgram(explodeProgram);
  delete [] _texture_data;
}



// Called to draw scene
void RenderScene(void)
{

//  static CStopWatch rotTimer;

  // Clear the window and the depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /** Set the modelView matrix */
  float viewDist = 5.0;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, -viewDist,  0, 0, 0,  0, 1, 0);
  glPushMatrix();
    glRotatef( 45.0f,0.0,1.0,0.0);
    glRotatef(-25.0f,1.0,0.0,0.0);

    /** Set the current 3D texture data from level set voxels */
    glActiveTextureARB(GL_TEXTURE0_ARB);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, szX, szY, szZ, 0,
                        GL_RGBA, GL_UNSIGNED_BYTE, _texture_data);


    /** draw some static stuff */
    chernobylPointGrid();

  /** done with modelview */
  glPopMatrix();
//  if(1) {
//  /** rotate the scene for exploding torus */
//  modelViewMatrix.Rotate(rotTimer.GetElapsedSeconds() * 10.0f, 0.0f, 1.0f, 0.0f);
//  modelViewMatrix.Rotate(rotTimer.GetElapsedSeconds() * 13.0f, 1.0f, 0.0f, 0.0f);

//  GLfloat vEyeLight[] = { -100.0f, 100.0f, 100.0f };
//  GLfloat vAmbientColor[] = { 0.1f, 0.1f, 0.1f, 1.0f };
//  GLfloat vDiffuseColor[] = { 0.1f, 1.0f, 0.1f, 1.0f };
//  GLfloat vSpecularColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };

//  glUseProgram(explodeProgram);
//  glUniformMatrix4fv(locMVP, 1, GL_FALSE, transformPipeline.GetModelViewProjectionMatrix());
//  glUniformMatrix4fv(locMV, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());
//  glUniformMatrix3fv(locNM, 1, GL_FALSE, transformPipeline.GetNormalMatrix());

//  float push_out = sinf(rotTimer.GetElapsedSeconds() * 3.0f) * 0.1f + 0.2f;

//  glUniform1f(locPushOut, push_out);

//  torusBatch.Draw();
//  }

 // modelViewMatrix.PopMatrix();

  glutSwapBuffers();
  glutPostRedisplay();
}


void ChangeSize(int w, int h)
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  GLdouble size = (GLdouble)((w >= h) ? w : h) / 2.0, aspect;
  if (w <= h) {
      aspect = (GLdouble)h/(GLdouble)w;
      glOrtho(-size, size, -size*aspect, size*aspect, -100000.0, 100000.0);
  }
  else {
      aspect = (GLdouble)w/(GLdouble)h;
      glOrtho(-size*aspect, size*aspect, -size, size, -100000.0, 100000.0);
  }
  glScaled(aspect, aspect, 1.0);
  glGetDoublev(GL_PROJECTION_MATRIX, _ortho_m);

  glLoadIdentity();
  gluPerspective(40.0/*AOV*/,
                 (GLdouble) w/h /* aspect ratio */,
                 0.1 /* zNear */,
                 20.0 /* zFar */ );
  glGetDoublev(GL_PROJECTION_MATRIX, _persp_m);
}

///////////////////////////////////////////////////////////////////////////////
// Main entry point for GLUT based programs
int main(int argc, char* argv[])
{

  //  gltSetWorkingDirectory(argv[0]);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STENCIL);
  glutInitWindowSize(800, 600);
  glutCreateWindow("Tex3D_Experimental");
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);

  GLenum err = glewInit();
  if (GLEW_OK != err) {
    fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
    return 1;
  }

  Setup_RenderingContext();
  glutMainLoop();
  Shutdown_RenderingContext();
  return 0;
}
