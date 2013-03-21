#version 150

precision highp float;

layout (triangles) in;
layout (triangle_strip, max_vertices = 12) out; // set max num outputs possible

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
  int num_verts=gl_in.length(); // should be 3
  for (n = 0; n < num_verts; n++)
  { /** If this shader works correctly,
                things should spin and flipflop color...    */
      color   =  vec4(1.0,0.0,0.0,0.25);
      color.x = (float(push_out>0))*color.x;
      gl_Position = mvpMatrix * vec4( gl_in[n].gl_Position.xyz,
                                      gl_in[n].gl_Position.w);
      EmitVertex();
  }
  EndPrimitive();
}

