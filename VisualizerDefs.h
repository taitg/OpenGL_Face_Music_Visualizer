
#ifndef VISUALIZERDEFS_H
#define VISUALIZERDEFS_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>
#include "glm/glm.hpp"
#include "kissfft/kiss_fft.h"

#define PI			3.14159265359
#define TWOPI		PI*2
#define HALFPI		PI/2
#define MSHPS		65535
#define MAXVERTS	MSHPS*4
#define MAXINDS		MAXVERTS*3
#define MAXSHAPES	32
#define N			2048
#define HALFN		N/2
#define DOUBLEN		N*2
#define SQRTHALFN	sqrt(HALFN)
#define TWOPIOVER	TWOPI / (N-1)
#define NUMTEX		1
#define TEX0		"texture.png"

#define PEAK_RISE 0.0001f

using namespace std;
using namespace glm;

// OpenGL window values
struct MyGLWindow {
	GLFWwindow *window;
	int w, h;
	double mousex, mousey;
	int initW = 1920, initH = 1080;
	int samples = 8;
};

// boolean flags
struct MyFlags {
	bool drawMode = false, clicking = false;
	bool wireframe = false, antialiasing = true;
	bool lighting = true, explode = true;
	bool simpleNormals = false, switchOuterNorm = true;
	bool peakRotateP = true, peakRotateT = true;
	bool peakFov = true, animTess = false;
	bool textureSky = true, colourSky = true;
};

// audio input parameters
struct MyAudioInput {
	unsigned int sampleRate = 44100;
	unsigned int numChannels = 2;
	unsigned int channel = 0;
	unsigned int device = 0;
};

// display parameters
struct MyDisplayValues {
	int scene = 0, mainShape = 1, shape = 1, colours = 2;
	int xscaling = 1, yscaling = 1, fftmode = 1;
	float yAdjust = 1.f;
	float volumeAdjust = 3.25f / SQRTHALFN;
	float xInc = 0.1f;
	int skySphereSplits = 12;
	int srevSplits = 24;
	float srevMax = TWOPI;
	float srevInc = PI * 2 / srevSplits;
	float tessInner = 1.f, tessOuter = 1.f;
};

// animation parameters
struct MyAnimation {
	float peak = 0.f, avgEnergy = 0.f;
	float lastFrameTime = 0.f;
	float targetFps = 66.6f, rateAdjust = 1.f;
	float upFactor = 0.004f, downfactor = 0.001f;
	float colourAdjust = 0.f, colourSpeed = 0.01f;
};

// uniform data
struct MyUniforms {
	mat4 viewMatrix, projMatrix;
	GLint viewL, projL;
	GLfloat camPoint[3];
	GLint camLocL;
	GLint simpleNormL;
	GLint tiL, toL;
	GLint peakL;
	GLint lightingL;
	GLint lightLocL;
	GLint intensityL, ambientL;
	GLint diffRatioL, phongL;
	GLint specColourL;
	GLint explodeL;
	GLint textureSkyL;
	GLint drawL;
};

// b-spline parameters
struct MyBSpline {
	int k = 3;
	double uinc = 0.005;
	double knots[HALFN];
	dvec2 c[HALFN];
};

// camera parameters
struct MyCamera {
	float camP = HALFPI, camT = HALFPI;
	float camR = 10.f;
	float inverted = 1.f;
	float scrollSpeed = 1.f, zoomSpeed = 0.15f;
	float minDistance = 1.f, maxDistance = 500.f;
	float fovDegrees = 38.f;
	float fov = fovDegrees * PI / 180.f;
	float zNear = .1f, zFar = 1000.f;
};

// lighting parameters
struct MyLight {
	GLfloat loc[3] = { 0.f, 0.f, 20.f };
	float intensity = 0.5f, ambient = 0.7f;
	float diffRatio = 1.0f;
	float phong = 90.0f;
	float specColour[3] = { 0.5f, 0.5f, 0.5f };
};

// shaders data
struct MyShader {
	// OpenGL names for vertex and fragment shaders, shader program
	GLuint  vertex;
	GLuint  TCS;
	GLuint  TES;
	GLuint	geometry;
	GLuint  fragment;
	GLuint  program;

	// initialize shader and program names to zero (OpenGL reserved value)
	MyShader() : vertex(0), fragment(0), TCS(0), TES(0), geometry(0), program(0)
	{}
};

// texture data
struct MyTexture {
	GLuint textureID;
	GLuint target;
	int width;
	int height;

	// initialize object names to zero (OpenGL reserved value)
	MyTexture() : textureID(0), target(0), width(0), height(0)
	{}
};

// geometry data
struct MyGeometry {
	// OpenGL names for array buffer objects, vertex array object
	GLuint  vertexBuffer;
	GLuint  textureBuffer;
	GLuint  colourBuffer;
	GLuint  elementBuffer;
	GLuint  vertexArray;
	GLsizei elementCount;

	// initialize object names to zero (OpenGL reserved value)
	MyGeometry() : vertexBuffer(0), colourBuffer(0), textureBuffer(0),
					elementBuffer(0), vertexArray(0), elementCount(0)
	{}
};

// data arrays and indices
struct MyDataValues {
	float xVals[N], previous[HALFN], points[HALFN][2];
	float pointShape[MAXVERTS][2], freeShape[MAXVERTS][2];
	float history[HALFN][HALFN][2];
	float shapePoints[HALFN];
	vector<dvec2> shapeControl;
	int scSelected = -1;
	float scRadius = 0.01f;
	int randoms[HALFN];
	GLfloat vertices[MAXVERTS][3];
	GLfloat colours[MAXVERTS][3];
	GLfloat texCoords[MAXVERTS][3];
	unsigned indices[MAXINDS];
	MyTexture texture;

	int numPoints = 0, numIndices = 0;
	int startVertex = 0, startIndex = 0;
	int numShapeVerts = 0;
	int numHistory = 0, maxHistory = 10;
	int randomIndex = 0;
};

// fast fourier transform data
struct MyFFT {
	kiss_fft_cfg cfg;
	kiss_fft_cpx in[N], out[N];
	float inBuffer[DOUBLEN];
};

#endif