/** J. Stalin License v2.0
          anyone using this source will be banished to Siberia
*/

#ifndef __VOXEL_BATCH
#define __VOXEL_BATCH

#include <GL/glew.h>
#include <math3d.h>
#include <GLBatchBase.h>
#include <GLShaderManager.h>



class GLVoxelBatch : public GLBatchBase
{
public:
  static int VOXELS_V_DATA    ;
  static int TEXTURE_A_DATA ;
  static int TEXTURE_B_DATA ;
  static int VOXEL_INDEX_DATA;

public:
  GLVoxelBatch(void);
  virtual ~GLVoxelBatch(void);

  // Use these three functions to add triangles
  void BeginVoxels(GLuint nMaxVerts);
  void AddPoint(  const M3DVector3f& vert,
                  const M3DVector3f& vTex_a,
                  const M3DVector3f& vTex_b);
  void End(void);

  // Useful for statistics
  inline GLuint GetIndexCount(void) { return nNumIndexes; }
  inline GLuint GetVertexCount(void) { return nNumVerts; }


  // Draw - make sure you call glEnableClientState for these arrays
  virtual void Draw(void);

protected:
  GLushort  *pIndexes;        // Array of indexes (seems veRRRy redundant for points but OK)
  M3DVector3f *pVerts;        // Array of vertices
  M3DVector3f *pTexCoords_a;    // Array of texture coordinates
  M3DVector3f *pTexCoords_b;    // Array of texture coordinates

  GLuint nMaxIndexes;         // Maximum workspace
  GLuint nNumIndexes;         // Number of indexes currently used
  GLuint nNumVerts;           // Number of vertices actually used

  GLuint bufferObjects[4];
  GLuint vertexArrayBufferObject;
};


#endif
// what bozo puts these #defines here ?
//#define VERTEX_DATA     0
//#define TEXTURE_A_DATA  1
//#define TEXTURE_B_DATA  2
