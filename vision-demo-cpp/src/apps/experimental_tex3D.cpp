
#include <iostream>
#include <math.h>
#include <GL/glew.h>
#include <GL/glu.h>
#include <GL/glext.h>           /* GLXContext */
#include <assert.h>
#include <stdio.h>   /** printf(%bla) */
#include <stdarg.h> /** varyadic funcs */
#include <GL/glut.h>
#include <X11/Intrinsic.h>    /* Display, Window */

namespace {
static GLuint  explodeProgram; // The geometry shader 'exploder' program

/** if these are set during rendering, you can use them correctly in the shader */
static GLint   locMVP;         // ModelViewProjection matrix uniform
static GLint   locMV;          // ModelView matrix uniform
static GLint   locNM;          // Normal matrix uniform
static GLint   locPushOut;     // How far to push the geomery out

static GLint   szX=32;
static GLint   szY=32;
static GLint   szZ=32;

double _persp_m[16], _ortho_m[16];	// projection matrices


//////////////////////////////////////////////////////////////////////////
// Load the shader from the source text
void gltLoadShaderSrc(const char *szShaderSrc, GLuint shader)
{
  GLchar *fsStringPtr[1];

  fsStringPtr[0] = (GLchar *)szShaderSrc;
  glShaderSource(shader, 1, (const GLchar **)fsStringPtr, NULL);
}


////////////////////////////////////////////////////////////////
// Load the shader from the specified file. Returns false if the
// shader could not be loaded
bool gltLoadShaderFile(const char *szFile, GLuint shader)
{
  GLint shaderLength = 0;
  FILE *fp;
  int MAX_SHADER_LENGTH = 4096;
  char* shaderText = new char [MAX_SHADER_LENGTH];

  // Open the shader file
  fp = fopen(szFile, "r");
  if(fp != NULL)
  {
    // See how long the file is
    while (fgetc(fp) != EOF)
      shaderLength++;

    // Allocate a block of memory to send in the shader
    assert(shaderLength < MAX_SHADER_LENGTH);   // make me bigger!
    if(shaderLength > MAX_SHADER_LENGTH)
    {
      fclose(fp);
      return false;
    }

    // Go back to beginning of file
    rewind(fp);

    // Read the whole file in
    if (shaderText != NULL)
      fread(shaderText, 1, shaderLength, fp);

    // Make sure it is null terminated and close the file
    shaderText[shaderLength] = '\0';
    fclose(fp);
  }
  else
    return false;
  std::cout <<"read shader text: " << std::endl << std::string(shaderText) << std::endl;
  // Load the string
  gltLoadShaderSrc((const char *)shaderText, shader);

  delete [] shaderText;
  return true;
}

/** Load a pair of shaders, compile, and link together. Specify the complete
  * source text for each shader. After the shader names, specify the number
  * of attributes, followed by the index and attribute name of each attribute */
GLuint gltLoadShaderTripletWithAttributes(const char *szVertexShader,
                                          const char *szGeometryShader,
                                          const char *szFragmentShader, ...)
{
  // Temporary Shader objects
  GLuint hVertexShader = 0;
  GLuint hGeometryShader = 0;
  GLuint hFragmentShader = 0;
  GLuint hReturn = 0;
  GLint testVal;
  int iArgCount;
  va_list attributeList;
  char *szNextArg = NULL;
  char infoLog[1024];

  // Create shader objects
  hVertexShader = glCreateShader(GL_VERTEX_SHADER);
  if (gltLoadShaderFile(szVertexShader, hVertexShader) == false)
    goto failed;
  glCompileShader(hVertexShader);
  glGetShaderiv(hVertexShader, GL_COMPILE_STATUS, &testVal);
  if  (testVal == GL_FALSE) {
    glGetShaderInfoLog(hVertexShader, 1024, NULL, infoLog);
    goto failed;
  }

  // Geometry shader is optional
  if (szGeometryShader) {
    hGeometryShader = glCreateShader(GL_GEOMETRY_SHADER);
    if(gltLoadShaderFile(szGeometryShader, hGeometryShader) == false)
      goto failed;
    glCompileShader(hGeometryShader);
    glGetShaderiv(hGeometryShader, GL_COMPILE_STATUS, &testVal);
    if  (testVal == GL_FALSE) {
      std::cout << "Failed to compile geometry shader!" << std::endl;
      glGetShaderInfoLog(hGeometryShader, 1024, NULL, infoLog);
      goto failed;
    }
  } else {
    std::cout << "OK, no geometry shader. " << std::endl;
  }

  // Fragment shader is optional (transform feedback only)
  if (szFragmentShader) {
    hFragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    if (gltLoadShaderFile(szFragmentShader, hFragmentShader) == false)
      goto failed;
    glCompileShader(hFragmentShader);
    glGetShaderiv(hFragmentShader, GL_COMPILE_STATUS, &testVal);
    if  (testVal == GL_FALSE) {
      glGetShaderInfoLog(hFragmentShader, 1024, NULL, infoLog);
      goto failed;
    }
  }

  std::cout << "done loading and compiling shaders!" << std::endl;

  // Create the final program object, and attach the shaders
  hReturn = glCreateProgram();
  glAttachShader(hReturn, hVertexShader);
  if (szGeometryShader)
    glAttachShader(hReturn, hGeometryShader);
  if (szFragmentShader)
    glAttachShader(hReturn, hFragmentShader);

  // Now, we need to bind the attribute names to their specific locations
  va_start(attributeList, szFragmentShader);

  // Iterate over this argument list
  iArgCount = va_arg(attributeList, int);	// Number of attributes
  for(int i = 0; i < iArgCount; i++)
  {
    int index = va_arg(attributeList, int);
    szNextArg = va_arg(attributeList, char*);
    glBindAttribLocation(hReturn, index, szNextArg);
  }
  va_end(attributeList);

  // Attempt to link
  glLinkProgram(hReturn);

  // These are no longer needed
  glDeleteShader(hVertexShader);
  glDeleteShader(hGeometryShader);
  glDeleteShader(hFragmentShader);

  // Make sure link worked too
  glGetProgramiv(hReturn, GL_LINK_STATUS, &testVal);
  if(testVal == GL_FALSE) {
    char infoLog[1024];
    glGetProgramInfoLog(hReturn, 1024, NULL, infoLog);
    fprintf(stderr, infoLog);
    goto failed;
  } else {
    std::cout << "Great, shader chain LINK success " << std::endl;
  }

  return hReturn; // DONE, Good

failed:
  std::cout<<"FAILED loading shaders!"<<std::endl;
  glDeleteProgram(hReturn);
  glDeleteShader(hFragmentShader);
  glDeleteShader(hGeometryShader);
  glDeleteShader(hVertexShader);

  return 0;
}



}

using namespace std;


void  chernobylPointGrid( ) {
  //  glEnable(GL_FRAGMENT_PROGRAM_ARB);
  //  glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, _prog[0]);
  if(0) {
    glBegin( GL_QUADS );
      glColor4f(1.0f,1.0f,1.0f,0.5f);
      glTexCoord3f(0.0,0.0,0.0); glVertex3f(-1.0,-1.0,0.0); glNormal3f(0.0f,0.0f,1.0f);
      glTexCoord3f(1.0,0.0,0.0); glVertex3f(+1.0,-1.0,0.0); glNormal3f(0.0f,0.0f,1.0f);
      glTexCoord3f(1.0,1.0,0.0); glVertex3f(+1.0,+1.0,0.0); glNormal3f(0.0f,0.0f,1.0f);
      glTexCoord3f(0.0,1.0,0.0); glVertex3f(-1.0,+1.0,0.0); glNormal3f(0.0f,0.0f,1.0f);
    glEnd();
  } else {
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
          glNormal3f(0.0f,0.0f,1.0f);
          glTexCoord3f( (x+1.0)/2.0, (y+1.0)/2.0, (z+1.0)/2.0 );
        }
      }
    }
    glEnd();
}

  // if drawing something else, disable texture 3D first!
  //  glDisable(GL_TEXTURE_3D);
  //  glDisable(GL_FRAGMENT_PROGRAM_ARB);


}

enum GLT_SHADER_ATTRIBUTE { GLT_ATTRIBUTE_VERTEX = 0,
                            GLT_ATTRIBUTE_COLOR,
                            GLT_ATTRIBUTE_NORMAL,
                            GLT_ATTRIBUTE_TEXTURE0, GLT_ATTRIBUTE_TEXTURE1,
                            GLT_ATTRIBUTE_TEXTURE2, GLT_ATTRIBUTE_TEXTURE3,
                            GLT_ATTRIBUTE_LAST};
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
        _texture_data[index+0]=255*(x<0); // R
        _texture_data[index+1]=255*(y>0); // G
        _texture_data[index+2]=255*(x>0); // B
        _texture_data[index+3]=128; // A
        index+=4;
        //(j*spectrum_width*tex_steps+i*spectrum_width+k)*4;
      }
    }
  }

  // Background
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f );
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

  explodeProgram = gltLoadShaderTripletWithAttributes("vertexX.shader",
                                                      NULL,//"geometryX.shader",
                                                      "fragmentX.shader",
                                                      2,
                                                      GLT_ATTRIBUTE_VERTEX, "vVertex",
                                                      GLT_ATTRIBUTE_NORMAL, "vNormal");

  locMVP = glGetUniformLocation(explodeProgram, "mvpMatrix");
  locMV  = glGetUniformLocation(explodeProgram, "mvMatrix");
  locNM  = glGetUniformLocation(explodeProgram, "normalMatrix");
  locPushOut = glGetUniformLocation(explodeProgram, "push_out");

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



/** lighting params, used in shaders (?) */
GLfloat vEyeLight[] = { -100.0f, 100.0f, 100.0f };
GLfloat vAmbientColor[] = { 0.1f, 0.1f, 0.1f, 1.0f };
GLfloat vDiffuseColor[] = { 0.1f, 1.0f, 0.1f, 1.0f };
GLfloat vSpecularColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
int steps = 0;
// Called to draw scene
void RenderScene(void)
{

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
  glEnable(GL_TEXTURE_3D);
  glActiveTextureARB(GL_TEXTURE0_ARB);  // the ARB stuff is probably not compatible?
  glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, szX, szY, szZ, 0,
               GL_RGBA, GL_UNSIGNED_BYTE, _texture_data);

  glDisable(GL_DEPTH_TEST);
//  glDisable(GL_TEXTURE_3D);
//  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE,GL_ONE);

  /** set shader pipeline to use now */
//  glUseProgram(explodeProgram);

//  /** set shader variables */
//  float push_out = sin(steps*0.01)*0.5  + 0.5;
//  glUniform1f(locPushOut, push_out);

//  float m[16];
//  glGetFloatv(GL_PROJECTION_MATRIX,m);
//  glUniformMatrix4fv(locMVP,1,GL_FALSE,m);

//  float m2[16];
//  glGetFloatv(GL_MODELVIEW_MATRIX,m2);
//  glUniformMatrix4fv(locMV,1,GL_FALSE,m2);


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


