#version 150

in vec3 position;
in vec2 textCoord;

out vec2 tc;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform int mode = 0;

void main()
{

    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0f);
    tc = textCoord;

}

