
// Vertex shader

#version 410

layout(location = 0) in vec3 VertexPosition;
layout(location = 1) in vec3 VertexColour;
layout(location = 2) in vec3 VertexTexture;

out vec3 vColour;
out vec3 vTexCoords;
out vec3 vPosition;

void main()
{
	// pass on data
	gl_Position = vec4(VertexPosition, 1.0);
    vColour = VertexColour;
	vTexCoords = VertexTexture;
}
