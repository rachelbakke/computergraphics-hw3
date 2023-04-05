#version 150

in vec3 position;
in vec3 normal;

out vec3 viewPosition;
out vec3 viewNormal;

uniform mat4 modelViewMatrix;
uniform mat4 normalMatrix;
uniform mat4 projectionMatrix;
uniform int mode = 0;

void main()
{
  // compute the transformed and projected vertex position (into gl_Position) 
  // compute the vertex color (into col)

  vec4 viewPosition4 = modelViewMatrix * vec4(normal, 1.0f);
  viewPosition = viewPosition4.xyz;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0f);
    viewNormal = normalize((normalMatrix * vec4(normal, 1.0f)).xyz);
}

