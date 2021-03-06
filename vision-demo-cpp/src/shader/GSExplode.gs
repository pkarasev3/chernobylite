// Geometry Shader 'Exploder' Example - Geometry Shader
// Graham Sellers
// OpenGL SuperBible 5th Edition
#version 150

precision highp float;

layout (triangles) in;
layout (triangle_strip, max_vertices = 12) out; // set max num outputs possible
/** Not allowed to mix output types!
      layout (points, max_vertices = 3) out; */

in Vertex
{
    vec3 normal;
    vec4 color;
} vertex[];

out vec4 color;

uniform vec3 vLightPosition;
uniform mat4 mvpMatrix;
uniform mat4 mvMatrix;
uniform mat3 normalMatrix;

uniform float push_out;

void main(void)
{
  int n;
  vec3 Fn =
  normalize(cross(gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz,
                  gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz));

  for (n = 0; n < gl_in.length(); n++) {
      color = vec4(1.0,0.0,0.0,0.25);//( n ) * vertex[n].color; //((float)n) / gl_in.length()
      gl_Position = mvpMatrix *
                     vec4( gl_in[n].gl_Position.xyz + Fn * push_out,
                                                 gl_in[n].gl_Position.w);
      EmitVertex();
  }
  EndPrimitive();

  // draw another triangle
  for (n = 0; n < gl_in.length(); n++) {
      color = vertex[n].color; //((float)n) / gl_in.length()
      gl_Position = mvpMatrix *
                     vec4( gl_in[n].gl_Position.xyz + Fn * (push_out+0.01f),
                                                 gl_in[n].gl_Position.w);
      EmitVertex();
  }
  EndPrimitive();


  // draw another triangle, pulsating at original position
  for (n = 0; n < gl_in.length(); n++) {
      color = vec4(0.0,1.0,0.1,0.25);//( n ) * vertex[n].color; //((float)n) / gl_in.length()
      //color = vertex[n].color; //((float)n) / gl_in.length()
      gl_Position = /** mvpMatrix *    For this one, show original position! */
                     vec4( gl_in[n].gl_Position.xyz + Fn * (push_out+0.01f),
                                                 gl_in[n].gl_Position.w);
      EmitVertex();
  }
  EndPrimitive();


}
