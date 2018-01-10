
// Tessellation control shader

#version 410
layout(vertices = 3) out;

in vec3 vColour[];
in vec3 vTexCoords[];

out vec3 tcColour[];
out vec3 tcTexCoords[];
out vec3 tcNormal[];

uniform float TessLevelInner;
uniform float TessLevelOuter;

#define ID gl_InvocationID

void main() {

	// set tessellation levels
    if (ID == 0) {
        gl_TessLevelInner[0] = TessLevelInner;
        gl_TessLevelOuter[0] = TessLevelOuter;
        gl_TessLevelOuter[1] = TessLevelOuter;
        gl_TessLevelOuter[2] = TessLevelOuter;
    }

	// pass on data
	gl_out[ID].gl_Position = gl_in[ID].gl_Position;		
    tcColour[ID] = vColour[ID]; 				
	tcTexCoords[ID] = vTexCoords[ID];

	// determine surface normal
	mat4 mod = mat4(1.0);
	vec4 c = vec4(0.0, 0.0, 0.0, 1.0);
	tcNormal[ID] = normalize(gl_out[ID].gl_Position.xyz - c.xyz);
}