// OpenGL toolkit
#include "GLTools.h"
#include "GLMatrixStack.h"
#include "GLFrame.h"
#include "GLFrustum.h"
#include "GLGeometryTransform.h"
#include "StopWatch.h"
#include "GLVoxelBatch.h"

#include <iostream>
#include <math.h>
#include <GL/glut.h>

using namespace std;

static GLFrame              viewFrame;
static GLFrustum            viewFrustum;
static GLMatrixStack        modelViewMatrix;
static GLMatrixStack        projectionMatrix;
static GLGeometryTransform  transformPipeline;
static GLShaderManager      shaderManager;

static GLVoxelBatch      voxelBatch;

static GLuint  explodeProgram; // The geometry shader 'exploder' program
static GLint   locMVP;         // ModelViewProjection matrix uniform
static GLint   locMV;          // ModelView matrix uniform
static GLint   locNM;          // Normal matrix uniform

static GLint   locPushOut;     // How far to push the geomery out

int szX=32; int szY=32; int szZ=32;

void  createChernobylPointGrid( ) {

  voxelBatch.BeginVoxels(szX*szY*szZ);
// argh , fixed-function shit, don't use
//  glDisable(GL_DEPTH_TEST);
//  glPointSize(3.0f);
//  glColor4f(1.0f,0.5f,0.5f,1.0f);
//  glBegin(GL_POINTS); //starts drawing of points
  for( int i=0; i<szY; i++ ) {
    float y=i*2.0/(szY-1)-1.0f;
    for( int j=0; j<szX; j++ ) {
      float x=j*2.0/(szX-1)-1.0f;
      for( int k=0; k<szZ; k++ ) {
        float z=k*2.0/(szZ-1)-1.0f;
        //glVertexAttrib1fv(GL_VERTEX_A)
//        glVertex3f(x,y,z);
//        glNormal3f(0.0f,0.0f,1.0f);
//        glTexCoord3f( (x+1.0)/2.0, (y+1.0)/2.0, (z+1.0)/2.0 );
        M3DVector3f verts{x*10,y*10,z*10};
        M3DVector3f vTex_a{0,0,0};
        M3DVector3f vTex_b{0,0,0};
        voxelBatch.AddPoint(&verts,&vTex_a,&vTex_b);
      }
    }
  } /** End the "chernobylGrid", wtf that is */
 // glEnd();
}

void drawChernobylPointGrid() {
  // TODO: call me !
   glPointSize(5.0);
   voxelBatch.Draw();

}

// This function does any needed initialization on the rendering context.
void SetupRC(void)
{
  // Background
  glClearColor(0.2f, 0.2f, 0.3f, 1.0f );

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_PROGRAM_POINT_SIZE);

  shaderManager.InitializeStockShaders();
  viewFrame.MoveForward(0.0f);

  // Make the torus
  //  gltMakeTorus(torusBatch, .70f /*R1*/, 0.10f/*r2*/, 11, 7);

  createChernobylPointGrid();

  explodeProgram = gltLoadShaderPair("vertexP.shader","fragmentX.shader");

//      gltLoadShaderTripletWithAttributes("vertexP.shader",
//                                                      NULL,// "geometryX.shader",
//                                                      "fragmentX.shader",
//                                                      1,
//                                                      GLT_ATTRIBUTE_VERTEX, "vVertex");

  locMVP = glGetUniformLocation(explodeProgram, "mvpMatrix");
  locMV  = glGetUniformLocation(explodeProgram, "mvMatrix");
  locNM  = glGetUniformLocation(explodeProgram, "normalMatrix");
  locPushOut = glGetUniformLocation(explodeProgram, "push_out");
}

// Cleanup
void ShutdownRC(void)
{
  glDeleteProgram(explodeProgram);
}

// Called to draw scene
void RenderScene(void)
{
  static CStopWatch rotTimer;

  // Clear the window and the depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//  viewFrame.SetOrigin(1.0,0.5,5.0);
//  viewFrame.SetForwardVector(-1,-0.5,-5);
//  viewFrame.SetUpVector(0,1,0);

//  viewFrame.MoveForward(-0.001f);
//  viewFrame.MoveRight(   0.0008f);
//   viewFrame.MoveUp(      0.01f);
//   viewFrame.MoveRight( 0.01f);
  //viewFrame.RotateLocal(0.5f,1.0,0.0,0.0);

  float xo = viewFrame.GetOriginX();
  float yo = viewFrame.GetOriginY();
  float zo = viewFrame.GetOriginZ();
//  cout << "xo= "<<xo << ",yo= "<<yo<<",zo= " << zo << "\r\n"; cout.flush();

  modelViewMatrix.LoadIdentity();

  modelViewMatrix.PushMatrix(viewFrame); // ahh, this does NOT change vertices going into GLSL

  modelViewMatrix.Rotate(rotTimer.GetElapsedSeconds() * 3.0f,0.0,1.0,0.0);

  // OK follow their style: manage your own matrices, don't use gl*** stack functions
//  glMatrixMode(GL_PROJECTION);
//  glLoadIdentity();
//  gluLookAt(0,0,5.0,0,0,0,0,1,0);       // but this does!

//  glMatrixMode(GL_MODELVIEW);
//  glLoadIdentity();



  if(1) {
    /** rotate the scene for exploding torus */
    // AHhhhhhh this stupid thing changes its internal matrix representation!
    //modelViewMatrix.Rotate(rotTimer.GetElapsedSeconds() * 10.0f, 0.0f, 1.0f, 0.0f);
    //modelViewMatrix.Rotate(rotTimer.GetElapsedSeconds() * 13.0f, 1.0f, 0.0f, 0.0f);

    GLfloat vEyeLight[] = { -100.0f, 100.0f, 100.0f };
    GLfloat vAmbientColor[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    GLfloat vDiffuseColor[] = { 0.1f, 1.0f, 0.1f, 1.0f };
    GLfloat vSpecularColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };

    glUseProgram(explodeProgram);

    glUniformMatrix4fv(locMVP, 1, GL_FALSE, transformPipeline.GetModelViewProjectionMatrix());

    glUniformMatrix4fv(locMV, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());

    glUniformMatrix3fv(locNM, 1, GL_FALSE, transformPipeline.GetNormalMatrix());

    float push_out = sinf(rotTimer.GetElapsedSeconds() * 3.0f) * 0.1f + 0.2f;

    glUniform1f(locPushOut, push_out);


    /** experimental, PK */
    drawChernobylPointGrid();

    //torusBatch.Draw();
  }

  //glPushMatrix();
//    glRotatef(+45.0,0.0,1.0,0.0);
//    glRotatef(-20.0,1.0,0.0,0.0);
    /** draw some static stuff */
  //glPopMatrix();
  modelViewMatrix.PopMatrix();  //   x>_x



  glutSwapBuffers();
  glutPostRedisplay();
}

void ChangeSize(int w, int h)
{
  // Prevent a divide by zero
  if(h == 0) h = 1;

  // Set Viewport to window dimensions
  glViewport(0, 0, w, h);

  viewFrustum.SetPerspective(60.0f, float(w)/float(h), 0.1f, 100.0f);

  projectionMatrix.LoadMatrix(viewFrustum.GetProjectionMatrix());
  transformPipeline.SetMatrixStacks(modelViewMatrix, projectionMatrix);
}

///////////////////////////////////////////////////////////////////////////////
// Main entry point for GLUT based programs
int main(int argc, char* argv[])
{
  gltSetWorkingDirectory(argv[0]);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STENCIL);
  glutInitWindowSize(800, 600);
  glutCreateWindow("Geometry Shader Exploder");
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);

  GLenum err = glewInit();
  if (GLEW_OK != err) {
    fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
    return 1;
  }

  SetupRC();
  glutMainLoop();
  ShutdownRC();
  return 0;
}
