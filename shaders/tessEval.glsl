
// Tessellation evaluation shader

#version 410
layout(triangles, equal_spacing, ccw) in;

in vec3 tcColour[];		
in vec3 tcTexCoords[];
in vec3 tcNormal[];

out vec3 teColour;		
out vec3 teTexCoords;
out vec3 tePosition;
out vec3 teNormal;

uniform mat4 view;
uniform mat4 proj;
uniform float peak;
uniform bool explode;
uniform bool draw;

void main()
{
	// determine vertex
    vec3 p0 = gl_TessCoord.x * gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_TessCoord.y * gl_in[1].gl_Position.xyz;
    vec3 p2 = gl_TessCoord.z * gl_in[2].gl_Position.xyz;
	vec3 newp = p0 + p1 + p2;

	// determine surface normal
	vec3 v1 = gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz;
	vec3 v2 = gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
	vec3 newn = normalize(cross(v2,v1));
	teNormal = newn;

	// move vertices
	if (tcTexCoords[0].z == 1) {
		float l0 = gl_TessCoord.x * length(gl_in[0].gl_Position.xy);
		float l1 = gl_TessCoord.y * length(gl_in[1].gl_Position.xy);
		float l2 = gl_TessCoord.z * length(gl_in[2].gl_Position.xy);
		float l = l0 + l1 + l2;
		vec3 newp2 = l * normalize(vec3(newp.x,newp.y,0));
		newp.x = newp2.x;
		newp.y = newp2.y;
	}
	else {
		float l0 = gl_TessCoord.x * length(gl_in[0].gl_Position.yz);
		float l1 = gl_TessCoord.y * length(gl_in[1].gl_Position.yz);
		float l2 = gl_TessCoord.z * length(gl_in[2].gl_Position.yz);
		float l = l0 + l1 + l2;
		vec3 newp2 = l * normalize(vec3(0,newp.y,newp.z));
		newp.y = newp2.y;
		newp.z = newp2.z;
		if (explode) newp += 0.2 * peak * l * newn;
	}

	// set position
	tePosition = newp;
	if (draw) 
		gl_Position = vec4(newp, 1.0);
	else
		gl_Position = proj * view * vec4(newp,1.0);
	//tePosition = gl_Position.xyz;

	// determine colour
	vec3 c0 = gl_TessCoord.x * tcColour[0];
	vec3 c1 = gl_TessCoord.y * tcColour[1];
	vec3 c2 = gl_TessCoord.z * tcColour[2];
	vec3 newc = c0 + c1 + c2;
	teColour = newc;

	// determine texCoords
	vec2 t0 = gl_TessCoord.x * tcTexCoords[0].xy;
	vec2 t1 = gl_TessCoord.y * tcTexCoords[1].xy;
	vec2 t2 = gl_TessCoord.z * tcTexCoords[2].xy;
	vec2 newtx = t0 + t1 + t2;
	teTexCoords = vec3(newtx,tcTexCoords[0].z);
}
