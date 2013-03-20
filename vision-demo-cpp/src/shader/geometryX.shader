precision highp float;

layout (points) in;
layout (triangle_strip, max_vertices = 12) out;

in Vertex
{
    vec3 normal;
    vec4 color;
} vertex[];

out vec4 color; /** gl_Position is implicitly an "out" as well*/

uniform vec3  vLightPosition;
uniform mat4  mvpMatrix;
uniform mat4  mvMatrix;
uniform mat3  normalMatrix;
uniform float push_out;

void main(void)
{
  int n;
  for (n = 0; n < 3; n++) {
      color       = vertex[0].color;

      gl_Position = gl_in[0].gl_Position+vec4(0.0,0.0,0.0,0.0);
      EmitVertex();

      gl_Position = gl_in[0].gl_Position+vec4(0.05,0.05,0.0,0.0);
      EmitVertex();

      gl_Position = gl_in[0].gl_Position+vec4(0.05,-0.05,0.05,0.0);
      EmitVertex();
  }
  EndPrimitive();

//          mvpMatrix *
//                     vec4( gl_in[n].gl_Position.xyz + Fn * (push_out+0.01f),
//                                                 gl_in[n].gl_Position.w);

}
