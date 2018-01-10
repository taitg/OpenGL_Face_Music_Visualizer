
// Fragment shader

#version 410

const float PI = 3.1415926535897932384626433832795;

in vec3 gColour;
in vec3 gTexCoords;
in vec3 tePosition;
in vec3 teNormal;
in vec3 gNormal;

out vec4 FragColour;

// textures
uniform sampler2D tex0;
uniform bool textureSky;

// light source
uniform bool showLight;
uniform vec3 light;
uniform vec3 specColour;
uniform float ambient;
uniform float diffRatio;
uniform float intensity;
uniform float phong;

// camera location
uniform vec3 camPoint;

// apply lighting model
vec4 applyLighting(vec4 colour) {

	vec4 newColour = colour;
	if (showLight) {
		vec3 norm = gNormal;
		vec3 pos = tePosition;

		vec3 camlight = camPoint;
		vec3 lightDir = -normalize(camlight - pos);
		vec3 viewRay = normalize(pos - camPoint);
		vec3 h = normalize(lightDir - viewRay);

		// ambient lighting
		newColour *= ambient;

		// diffuse lighting
		newColour += colour * diffRatio * intensity * max(0.0, dot(norm, lightDir));

		// specular lighting
		float maxTerm = max(0.0, dot(norm, h));
		vec3 specular = specColour * intensity * pow(maxTerm, phong);
		newColour += vec4(specular, 1.0);
	}
	return newColour;
}

// main function
void main(void) {

	vec4 colour;
	
	if (textureSky && gTexCoords.z == 1) 
		colour = 0.25*colour + 0.75*applyLighting(texture(tex0, gTexCoords.xy));
	else 
		colour = applyLighting(vec4(gColour,1.0));
	
	FragColour = colour;
}
