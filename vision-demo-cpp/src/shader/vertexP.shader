// Geometry Shader 'Exploder' Example - Vertex Shader
// Graham Sellers
// OpenGL SuperBible 5th Edition
#version 150

precision highp float;

// Incoming per vertex... position 
in vec4 vVertex;
//in vec3 vNormal;

//out Vertex
/*{
    vec3 normal;
    vec4 color;
} vertex;
*/

out color;

uniform vec3 vLightPosition;
uniform mat4 mvMatrix;
uniform mat3 normalMatrix;

void main(void)
{
    vec3 vNormal = vec3(0,0,1);
    // Get surface normal in eye coordinates
    vec3 vEyeNormal = normalMatrix * vNormal;

    // Get vertex position in eye coordinates
    vec4 vPosition4 = mvMatrix * vVertex; // mvMatrix is bullshit here
    vec3 vPosition3 = vPosition4.xyz / vPosition4.w;

    // Get vector to light source
    vec3 vLightDir = normalize(vLightPosition - vPosition3);

    gl_PointSize = 2.0;

    // Dot product gives us diffuse intensity
    color = vec4( 0.0, 
                         1.0, 
                         1.0, 
                         1.0);
// mvMatrix * 
    gl_Position   = mvMatrix * vVertex; //mvMatrix * 
    //vertex.normal = vNormal;
}
