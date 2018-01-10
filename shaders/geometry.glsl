
// Geometry shader

#version 410
layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in vec3 teColour[3];
in vec3 teTexCoords[3];
in vec3 tePosition[3];
in vec3 teNormal[3];

out vec3 gColour;
out vec3 gTexCoords;
out vec3 gNormal;

uniform bool simpleNorm;

void main() {

	// calculate normal
	vec3 A = tePosition[2] - tePosition[0];
	vec3 B = tePosition[1] - tePosition[0];
	vec3 norm = normalize(cross(B, A));
    
	// vertex 1
	gColour = teColour[0];
	gTexCoords = teTexCoords[0];
    gl_Position = gl_in[0].gl_Position;
	if (simpleNorm) gNormal = teNormal[0];
	else gNormal = norm;
	EmitVertex();

	// vertex 2
	gColour = teColour[1];
	gTexCoords = teTexCoords[1];
    gl_Position = gl_in[1].gl_Position; 
	if (simpleNorm) gNormal = teNormal[1];
	else gNormal = norm;
	EmitVertex();

	// vertex 3
	gColour = teColour[2];
	gTexCoords = teTexCoords[2];
    gl_Position = gl_in[2].gl_Position; 
	if (simpleNorm) gNormal = teNormal[2];
	else gNormal = norm;
	EmitVertex();

    EndPrimitive();
}