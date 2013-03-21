// Geometry Shader 'Explode' Example
// OpenGL SuperBible 5th Edition
// Demonstrates 'explosion' of geometry using OpenGL geometry shaders
// Program by Graham Sellers

// OpenGL toolkit
#include "GLTools.h"
#include "GLMatrixStack.h"
#include "GLFrame.h"
#include "GLFrustum.h"
#include "GLGeometryTransform.h"
#include "StopWatch.h"
#include <iostream>

#include <math.h>
#ifdef __APPLE__
#include <glut/glut.h>
#else
#define FREEGLUT_STATIC
#include <GL/glut.h>
#endif

#ifdef _MSC_VER
#pragma comment (lib, "GLTools.lib")
#endif /* _MSC_VER */
#include "GLVoxelBatch.h"
using namespace std;

static GLFrame              viewFrame;
static GLFrustum            viewFrustum;
static GLTriangleBatch      torusBatch;
static GLMatrixStack        modelViewMatrix;
static GLMatrixStack        projectionMatrix;
static GLGeometryTransform  transformPipeline;
static GLShaderManager      shaderManager;

static GLuint  explodeProgram; // The geometry shader 'exploder' program
static GLint   locMVP;         // ModelViewProjection matrix uniform
static GLint   locMV;          // ModelView matrix uniform
static GLint   locNM;          // Normal matrix uniform

static GLint   locPushOut;     // How far to push the geomery out

static GLVoxelBatch  voxelBatch;

int szX=16; int szY = 16; int szZ=16;

void makeVoxelGrid() {
  voxelBatch.BeginVoxels(szX*szY*szZ);

  for( int i=0; i<szY; i++ ) {
    float y=i*2.0/(szY-1)-1.0f;
    for( int j=0; j<szX; j++ ) {
      float x=j*2.0/(szX-1)-1.0f;
      for( int k=0; k<szZ; k++ ) {
        float z=k*2.0/(szZ-1)-1.0f;
        M3DVector3f verts{x,y,z};
        M3DVector3f vTex_a{0,0,0};
        M3DVector3f vTex_b{0,0,0};
        voxelBatch.AddPoint(verts,vTex_a,vTex_b);
      }
    }
  }
}

// This function does any needed initialization on the rendering context.
void SetupRC(void)
{
  // Background
  glClearColor(0.2f, 0.2f, 0.3f, 1.0f );

  glEnable(GL_DEPTH_TEST);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_PROGRAM_POINT_SIZE);
  shaderManager.InitializeStockShaders();
  viewFrame.MoveForward(4.0f);

  // Make the torus
  gltMakeTorus(torusBatch, .70f /*R1*/, 0.10f/*r2*/, 11, 7);

  makeVoxelGrid();

  explodeProgram = gltLoadShaderTripletWithAttributes("vertexP.shader",
                                                      "geometryX.shader",
                                                      "fragmentX.shader",
                                                          2,
                                                          GLT_ATTRIBUTE_VERTEX, "vVertex",
                                                          GLT_ATTRIBUTE_NORMAL, "vNormal");

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

  modelViewMatrix.PushMatrix(viewFrame);

  if(1) {
    /** rotate the scene for exploding torus */
    modelViewMatrix.Rotate(rotTimer.GetElapsedSeconds() * 10.0f, 0.0f, 1.0f, 0.0f);
    modelViewMatrix.Rotate(rotTimer.GetElapsedSeconds() * 13.0f, 1.0f, 0.0f, 0.0f);

    GLfloat vEyeLight[] = { -100.0f, 100.0f, 100.0f };
    GLfloat vAmbientColor[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    GLfloat vDiffuseColor[] = { 0.1f, 1.0f, 0.1f, 1.0f };
    GLfloat vSpecularColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };

    glUseProgram(explodeProgram);

    glUniformMatrix4fv(locMVP, 1, GL_FALSE, transformPipeline.GetModelViewProjectionMatrix());
    //  float m[16];
    //  glGetFloatv(GL_PROJECTION_MATRIX,m);
    //  glUniformMatrix4fv(locMVP,1,GL_FALSE,m);

    //glUniformMatrix4fv(locMV, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());
    float m2[16];
    glGetFloatv(GL_MODELVIEW_MATRIX,m2);
    glUniformMatrix4fv(locMV,1,GL_FALSE,m2);

    glUniformMatrix3fv(locNM, 1, GL_FALSE, transformPipeline.GetNormalMatrix());

    float push_out = sinf(rotTimer.GetElapsedSeconds() * 20.0f);
    glUniform1f(locPushOut, push_out);

    torusBatch.Draw();
    //voxelBatch.Draw();
  }
  modelViewMatrix.PopMatrix();

  torusBatch.Draw();

  glutSwapBuffers();
  glutPostRedisplay();
}

void ChangeSize(int w, int h)
{
  // Prevent a divide by zero
  if(h == 0) h = 1;

  // Set Viewport to window dimensions
  glViewport(0, 0, w, h);

  viewFrustum.SetPerspective(35.0f, float(w)/float(h), 1.0f, 100.0f);

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
  glutCreateWindow("Geometry Shader TEST");
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
