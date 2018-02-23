#version 330 core
// When you edit these shaders, Clear CMake Configuration so they are copied to build folders
// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/

in vec2 UV;
out vec3 FragColor;

uniform int hasNormals;
uniform int hasTextures;
uniform sampler2D myTextureSampler;

void main() {
    FragColor = texture(myTextureSampler, UV).rgb;
}
