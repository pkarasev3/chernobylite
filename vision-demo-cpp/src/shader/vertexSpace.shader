// SpaceFlight Shader
// Vertex Shader
// Richard S. Wright Jr.
// OpenGL SuperBible
#version 130

// Incoming per vertex... position and normal
in vec4 vVertex;
in vec4 vColor;

uniform mat4   mvpMatrix;
uniform float  timeStamp;

out vec4 vStarColor;

void main(void)
    {
    vec4 vNewVertex = vVertex;
    vStarColor = vColor;

    // Offset by running time, makes it move closer
    //vNewVertex.z += timeStamp;
    vStarColor.b = max(0.0,min(1.0,(vNewVertex.z)));
    //vStarColor.r = max(0.0,min(1.0,(vNewVertex.x)));
    // If out of range, adjust
    //if(vNewVertex.z > -1.0)
    //    vNewVertex.z -= 999.0;

    gl_PointSize = 2.0;// + (vNewVertex.z / sqrt(-vNewVertex.z));

    // If they are very small, fade them up
//    if(gl_PointSize < 4.0)
//        vStarColor = smoothstep(0.0, 4.0, gl_PointSize) * vStarColor;


    // Don't forget to transform the geometry!
    gl_Position = mvpMatrix * vNewVertex;
    }
