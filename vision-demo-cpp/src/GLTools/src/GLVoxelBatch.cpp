/** J. Stalin License v2.0
          anyone using this source will be banished to Siberia
*/

#include <GLVoxelBatch.h>
#include <GLShaderManager.h>

int GLVoxelBatch::VOXELS_V_DATA  = 0;
int GLVoxelBatch::TEXTURE_A_DATA = 1;
int GLVoxelBatch::TEXTURE_B_DATA = 2;
int GLVoxelBatch::VOXEL_INDEX_DATA = 3;


///////////////////////////////////////////////////////////
GLVoxelBatch::GLVoxelBatch(void)
{

  pVerts       = NULL;
  pTexCoords_a = NULL;
  pTexCoords_b = NULL;

  nMaxIndexes = 0;
  nNumIndexes = 0;
  nNumVerts   = 0;
}

////////////////////////////////////////////////////////////
// Free any dynamically allocated memory. For those C programmers
// coming to C++, it is perfectly valid to delete a NULL pointer.
GLVoxelBatch::~GLVoxelBatch(void)
{
  // Just in case these still are allocated when the object is destroyed
  //delete [] pIndexes;
  delete [] pVerts;
  //delete [] pNorms;
  delete [] pTexCoords_a;
  delete [] pTexCoords_b;

  // Delete buffer objects
  glDeleteBuffers(4, bufferObjects);

#ifndef OPENGL_ES
  glDeleteVertexArrays(1, &vertexArrayBufferObject);
#endif
}

////////////////////////////////////////////////////////////
// Start assembling a mesh. You need to specify a maximum amount
// of indexes that you expect.
void GLVoxelBatch::BeginVoxels(GLuint nMaxVerts)
{
  // Just in case this gets called more than once...
  if( pVerts )
    delete [] pVerts;
  if( pTexCoords_a )
    delete [] pTexCoords_a;
  if( pTexCoords_b)
    delete [] pTexCoords_b;
  if( pIndexes )
    delete [] pIndexes;

  nMaxIndexes = nMaxVerts;
  nNumIndexes = 0;
  nNumVerts   = 0;

  // Allocate new blocks. In reality, the other arrays will be
  // much shorter than the index array
  pIndexes     = new GLushort[nMaxIndexes];
  pVerts       = new M3DVector3f[nMaxIndexes];
  pTexCoords_a = new M3DVector3f[nMaxIndexes];
  pTexCoords_b = new M3DVector3f[nMaxIndexes];
  //    pNorms = new M3DVector3f[nMaxIndexes];
  //    pTexCoords = new M3DVector2f[nMaxIndexes];
}

/////////////////////////////////////////////////////////////////
// add another 3D voxel-grid-ish point
void GLVoxelBatch::AddPoint( const M3DVector3f& vert,
                             const M3DVector3f& vTex_a,
                             const M3DVector3f& vTex_b)
{
  const float _err = 0.00001f;

  //memcpy(pNorms[nNumVerts], vNorms[iVertex], sizeof(M3DVector3f));
  GLuint iVertex = 0; //iVertex < 3; iVertex++)
  memcpy(&(pVerts[nNumVerts]), &(vert[0]), sizeof(M3DVector3f));

  float tx=pVerts[nNumVerts][0];
  float ty=pVerts[nNumVerts][1];
  float tz=pVerts[nNumVerts][2];

  memcpy(pTexCoords_a[nNumVerts], &(vTex_a[0]), sizeof(M3DVector3f));
  memcpy(pTexCoords_b[nNumVerts], &(vTex_b[0]), sizeof(M3DVector3f));
  pIndexes[nNumIndexes] = nNumVerts;
  nNumIndexes++;
  nNumVerts++;
}



//////////////////////////////////////////////////////////////////
// Compact the data. This is a nice utility, but you should really
// save the results of the indexing for future use if the model data
// is static (doesn't change).
void GLVoxelBatch::End(void)
{
#ifndef OPENGL_ES
  // Create the master vertex array object
  glGenVertexArrays(1, &vertexArrayBufferObject);
  glBindVertexArray(vertexArrayBufferObject);
#endif

  // Create the buffer objects
  glGenBuffers(4, bufferObjects);

  // Copy data to video memory
  // Vertex data
  glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[VOXELS_V_DATA]);
  glEnableVertexAttribArray(GLT_ATTRIBUTE_VERTEX);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nNumVerts*3, pVerts, GL_STATIC_DRAW);
  glVertexAttribPointer(GLT_ATTRIBUTE_VERTEX, 3, GL_FLOAT, GL_FALSE, 0, 0);


//  // Normal data
//  glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[NORMAL_DATA]);
//  glEnableVertexAttribArray(GLT_ATTRIBUTE_NORMAL);
//  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nNumVerts*3, pNorms, GL_STATIC_DRAW);
//  glVertexAttribPointer(GLT_ATTRIBUTE_NORMAL, 3, GL_FLOAT, GL_FALSE, 0, 0);

  // Texture coordinates
  glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[TEXTURE_A_DATA]);
  glEnableVertexAttribArray(GLT_ATTRIBUTE_TEXTURE0);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nNumVerts*3, pTexCoords_a, GL_STATIC_DRAW);
  glVertexAttribPointer(GLT_ATTRIBUTE_TEXTURE0, 2, GL_FLOAT, GL_FALSE, 0, 0);

  // Texture coordinates
  glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[TEXTURE_B_DATA]);
  glEnableVertexAttribArray(GLT_ATTRIBUTE_TEXTURE1);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nNumVerts*2, pTexCoords_b, GL_STATIC_DRAW);
  glVertexAttribPointer(GLT_ATTRIBUTE_TEXTURE1, 2, GL_FLOAT, GL_FALSE, 0, 0);

  // Indexes
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferObjects[VOXEL_INDEX_DATA]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLushort)*nNumIndexes, pIndexes, GL_STATIC_DRAW);


  // Done
#ifndef OPENGL_ES
  glBindVertexArray(0);
#endif

  // Free older, larger arrays
  delete [] pIndexes;
  delete [] pVerts;
  delete [] pTexCoords_a;
  delete [] pTexCoords_b;

  // Reasign pointers so they are marked as unused
  pIndexes = NULL;
  pVerts = NULL;
  pTexCoords_a = NULL;
  pTexCoords_b = NULL;

  // Unbind to anybody
#ifndef OPENGL_ES
  glBindVertexArray(0);
#endif
}

//////////////////////////////////////////////////////////////////////////
// Draw - make sure you call glEnableClientState for these arrays
void GLVoxelBatch::Draw(void)
{
#ifndef OPENGL_ES
  glBindVertexArray(vertexArrayBufferObject);
#else
  glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[VERTEX_DATA]);
  glEnableVertexAttribArray(GLT_ATTRIBUTE_VERTEX);
  glVertexAttribPointer(GLT_ATTRIBUTE_VERTEX, 3, GL_FLOAT, GL_FALSE, 0, 0);

  // Normal data
  glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[NORMAL_DATA]);
  glEnableVertexAttribArray(GLT_ATTRIBUTE_NORMAL);
  glVertexAttribPointer(GLT_ATTRIBUTE_NORMAL, 3, GL_FLOAT, GL_FALSE, 0, 0);

  // Texture coordinates
  glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[TEXTURE_DATA]);
  glEnableVertexAttribArray(GLT_ATTRIBUTE_TEXTURE0);
  glVertexAttribPointer(GLT_ATTRIBUTE_TEXTURE0, 2, GL_FLOAT, GL_FALSE, 0, 0);

  // Indexes
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferObjects[INDEX_DATA]);
#endif


  glDrawElements(GL_POINTS, nNumIndexes, GL_UNSIGNED_SHORT, pIndexes);

#ifndef OPENGL_ES
  // Unbind to anybody
  glBindVertexArray(0);
#else
  glDisableVertexAttribArray(GLT_ATTRIBUTE_VERTEX);
  glDisableVertexAttribArray(GLT_ATTRIBUTE_NORMAL);
  glDisableVertexAttribArray(GLT_ATTRIBUTE_TEXTURE0);
#endif
}


