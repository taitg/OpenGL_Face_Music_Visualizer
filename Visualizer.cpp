
#include <Windows.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "glm/gtc/type_ptr.hpp"

#include <fstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cstdlib>
#include <string>
#include <math.h>
#include <ctime>
#include <random>
#include <map>
#include <vector>

#include "rtaudio/RtAudio.h"
#include "kissfft/kiss_fft.h"
#include "VisualizerDefs.h"
#include "GLsupport.h"

using namespace std;
using namespace glm;

// --------------------------------------------------------------------------
// Global vars (structs defined in VisualizerDefs.h)

MyGLWindow glWindow;
MyAudioInput ain;
MyFlags flags;
MyDataValues dataVals;
MyFFT fft;
MyDisplayValues display;
MyAnimation anim;
MyCamera camera;
MyBSpline bspline;
MyLight light;
HANDLE ghMutex;

// --------------------------------------------------------------------------
// Setter functions

/*
* Set the OpenGL window values
*/
void SetGLWindow(int initW, int initH, int samples) {

	glWindow.initW = initW;
	glWindow.initH = initH;
	glWindow.samples = samples;
}

/*
* Set the audio input values
*/
void SetAudioInput(int sampleRate, int numChannels, int channel, int device) {

	ain.sampleRate = sampleRate;
	ain.numChannels = numChannels;
	ain.channel = channel;
	ain.device = device;
}

/*
* Set the flags
*/
void SetFlags(bool draw, bool animTess, bool aa, bool explode, bool lighting,
				bool peakFov, bool peakRotate, bool simpleNorms,
				bool switchNorms, bool texSky, bool cSky, bool wireframe) {

	flags.drawMode = draw;
	flags.animTess = animTess;
	flags.antialiasing = aa;
	flags.explode = explode;
	flags.lighting = lighting;
	flags.peakFov = peakFov;
	flags.peakRotateP = peakRotate;
	flags.peakRotateP = peakRotate;
	flags.simpleNormals = simpleNorms;
	flags.switchOuterNorm = switchNorms;
	flags.textureSky = texSky;
	flags.colourSky = cSky;
	flags.wireframe = wireframe;
}

/*
* Set the display values
*/
void SetDisplayValues(float volumeAdjust, float xInc, int srevSplits,
						int ssSplits, float ti, float to,
						int colour, int shape) {

	display.volumeAdjust = volumeAdjust / SQRTHALFN;
	display.xInc = xInc;
	display.srevSplits = srevSplits;
	display.srevInc = PI * 2 / srevSplits;
	display.skySphereSplits = ssSplits;
	display.tessInner = ti;
	display.tessOuter = to;
	display.colours = colour;
	display.shape = shape;
}

/*
* Set the bspline values
*/
void SetBSpline(int k, float uinc) {

	bspline.k = k;
	bspline.uinc = uinc;
}

/*
* Set the camera values
*/
void SetCamera(float scrollSpeed, float zoomSpeed, float fovDegrees) {

	camera.scrollSpeed = scrollSpeed;
	camera.zoomSpeed = zoomSpeed;
	camera.fovDegrees = fovDegrees;
	camera.fov = fovDegrees * PI / 180.f;
}

/*
* Set lighting values
*/
void SetLight(float intensity, float ambient, float diffRatio, 
				float phong, float specColour[3]) {

	light.intensity = intensity;
	light.ambient = ambient;
	light.diffRatio = diffRatio;
	light.phong = phong;
	light.specColour[0] = specColour[0];
	light.specColour[1] = specColour[1];
	light.specColour[2] = specColour[2];
}

// --------------------------------------------------------------------------
/*
* Called when audio buffer is filled with data and ready,
* copies data to fft buffer
*/
int AudioCallback(void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
	double streamTime, RtAudioStreamStatus status, void *userData) {

	if (status) cout << "Stream overflow detected!" << endl;

	// Request ownership of mutex
	DWORD dwCount = 0, dwWaitResult;
	while (dwCount < 20) {
		dwWaitResult = WaitForSingleObject(ghMutex, INFINITE);  // no time-out interval
		switch (dwWaitResult) {
			// The thread got ownership of the mutex
			case WAIT_OBJECT_0:
				__try {

					// copy values to fft in-buffer
					memcpy(fft.inBuffer, inputBuffer, nBufferFrames * 2 * sizeof(float));
					for (int i = 0; i < N; i++) {

						// apply windowing function, fill in array
						float window = 0.5f * (1 - cos(i * TWOPIOVER));
						fft.in[i].r = window * fft.inBuffer[i];
						fft.in[i].i = 0;

						// set peak value
						if (fft.in[i].r > anim.peak) anim.peak += PEAK_RISE;//= fft.in[i].r;
					}
					dwCount++;
				}
				__finally {
					// Release ownership of the mutex object
					if (!ReleaseMutex(ghMutex)) {
						cout << "ERROR: could not release mutex" << endl;
						cin.get();
						exit(EXIT_FAILURE);
					}
				}
				break;
			// The thread got ownership of an abandoned mutex
			// The state is indeterminate
			case WAIT_ABANDONED: return 0;
		}
	}
	return 0;
}
// --------------------------------------------------------------------------
/*
* Generate a 2D point on the B-spline curve for a given u
*/
dvec2 GetBspline(float points[][2], double u, int d) {

	for (int i = 0; i <= bspline.k - 1; i++) {
		bspline.c[i].x = points[d - i][0];
		bspline.c[i].y = points[d - i][1];
	}

	for (int r = bspline.k; r >= 2; r--) {
		int i = d;
		for (int s = 0; s <= r - 2; s++) {
			double u_i = bspline.knots[i];
			double u_ir1 = bspline.knots[i + r - 1];
			double omega = (u - u_i) / (u_ir1 - u_i);
			bspline.c[s] = omega * bspline.c[s] + (1 - omega) * bspline.c[s + 1];
			i--;
		}
	}
	return bspline.c[0];
}

/*
* Generate a 2D point on the B-spline curve from a given u
* Takes a vector of control points instead of an array
*/
dvec2 GetBSplineVec(vector<dvec2> control, double u, int d) {

	dvec2 *cv = new dvec2[control.size()];
	for (int i = 0; i <= bspline.k - 1; i++) {
		cv[i] = control[d - i];
	}

	for (int r = bspline.k; r >= 2; r--) {
		int i = d;
		for (int s = 0; s <= r - 2; s++) {
			double u_i = bspline.knots[i];
			double u_ir1 = bspline.knots[i + r - 1];
			double omega = (u - u_i) / (u_ir1 - u_i);
			cv[s] = omega * cv[s] + (1 - omega) * cv[s + 1];
			i--;
		}
	}

	dvec2 result = cv[0];
	delete cv;
	return result;
}

/*
* Generate delta knot index for a given u
*/
int Delta(int numPoints, double u) {
	int max = numPoints + bspline.k - 1;
	for (int i = 0; i <= max; i++) {
		if (u >= bspline.knots[i] && u < bspline.knots[i + 1])
			return i;
	}
	return -1;
}

/*
* Generate a standard knot sequence based
* on the current order and number of control points
*/
void GenerateKnots(int numPoints) {

	for (int i = 0; i < bspline.k; i++)
		bspline.knots[i] = 0;

	int middle = numPoints - bspline.k;
	for (int i = 0; i < middle; i++)
		bspline.knots[i + bspline.k] = double(i + 1) / (middle + 1);

	for (int i = middle + bspline.k; i < middle + 2 * bspline.k; i++)
		bspline.knots[i] = 1;
}
// --------------------------------------------------------------------------
/*
* Get a random int in the range
*/
int RandomGen(int min, int max) {
	mt19937 rng;
	rng.seed(random_device()());
	uniform_int_distribution<mt19937::result_type> dist6(min, max);
	return dist6(rng);
}

/*
* Fill random array with ints
*/
void InitRandoms(int low1, int high1, int low2, int high2) {
	for (int i = 0; i + 1 < HALFN; i += 2) {
		dataVals.randoms[i] = RandomGen(low1, high1);
		dataVals.randoms[i+1] = RandomGen(low2, high2);
	}
	dataVals.randoms[HALFN] = low2;
}

/*
* Create a sky sphere using triangles
*/
void GenerateSkySphere() {

	float R = 100.f;
	float unit = PI / display.skySphereSplits;
	float x = 0.0, y = 0.0;
	int i = 0, j = 0;
	float px, py, pz;

	// iterate over rings
	for (float it = 0.f; it < PI; it += unit) {
		x = 0.0;

		// iterate over points of each ring
		for (float ip = 0.f; ip < TWOPI; ip += unit) {
			float nextp = ip + unit;
			float nextt = it + unit;

			// p1
			px = R * cos(ip) * sin(it);
			py = R * sin(ip) * sin(it);
			pz = R * cos(it);
			dataVals.vertices[i][0] = px;
			dataVals.vertices[i][1] = py;
			dataVals.vertices[i][2] = pz;
			dataVals.colours[i][0] = 0;
			dataVals.colours[i][1] = 0;
			dataVals.colours[i][2] = 0;
			dataVals.texCoords[i][0] = 1.0 - x;
			dataVals.texCoords[i][1] = y;
			dataVals.texCoords[i][2] = 1; // 1 = sky texture
			x += 1.0 / (2.0 * display.skySphereSplits);
			i++;

			// p2
			px = R * cos(nextp) * sin(it);
			py = R * sin(nextp) * sin(it);
			pz = R * cos(it);
			dataVals.vertices[i][0] = px;
			dataVals.vertices[i][1] = py;
			dataVals.vertices[i][2] = pz;
			dataVals.colours[i][0] = 0;
			dataVals.colours[i][1] = 0;
			dataVals.colours[i][2] = 0;
			dataVals.texCoords[i][0] = 1.0 - x;
			dataVals.texCoords[i][1] = y;
			dataVals.texCoords[i][2] = 1;
			i++;

			// p3
			px = R * cos(ip) * sin(nextt);
			py = R * sin(ip) * sin(nextt);
			pz = R * cos(nextt);
			dataVals.vertices[i][0] = px;
			dataVals.vertices[i][1] = py;
			dataVals.vertices[i][2] = pz;
			dataVals.colours[i][0] = 0;
			dataVals.colours[i][1] = 0;
			dataVals.colours[i][2] = 0;
			dataVals.texCoords[i][0] = 1.0 -
				(x - 1.0 / (2.0 * display.skySphereSplits));
			dataVals.texCoords[i][1] = y + 1.0 / display.skySphereSplits;
			dataVals.texCoords[i][2] = 1;
			i++;

			int offset = 0;
			if (ip > TWOPI - unit) {
				// p4
				px = R * cos(nextp) * sin(nextt);
				py = R * sin(nextp) * sin(nextt);
				pz = R * cos(nextt);
				dataVals.vertices[i][0] = px;
				dataVals.vertices[i][1] = py;
				dataVals.vertices[i][2] = pz;
				dataVals.colours[i][0] = 1;
				dataVals.colours[i][1] = 1;
				dataVals.colours[i][2] = 1;
				dataVals.texCoords[i][0] = 1.0 - x;
				dataVals.texCoords[i][1] = y + 1.0 / display.skySphereSplits;
				dataVals.texCoords[i][2] = 1;
				i++;
				offset = 1;
			}

			// indices
			dataVals.indices[j] = i - 1 - offset;
			j++;
			dataVals.indices[j] = i - 2 - offset;
			j++;
			dataVals.indices[j] = i - 3 - offset;
			j++;

			if (it >= 0.0 && it < PI - unit) {
				dataVals.indices[j] = i - 1 - offset;
				j++;
				if (flags.switchOuterNorm) {
					dataVals.indices[j] = i - 2 - offset;
					j++;
				}
				if (ip < TWOPI - unit) {
					dataVals.indices[j] = i + 2 - offset;
					j++;
				}
				else {
					dataVals.indices[j] = i - offset;
					j++;
				}
				if (!flags.switchOuterNorm) {
					dataVals.indices[j] = i - 2 - offset;
					j++;
				}
			}
		}
		y += 1.0 / display.skySphereSplits;
	}

	dataVals.startVertex = i;
	dataVals.startIndex = j;
}

/*
* Gets a sequence of control points from fft to use for bspline
*/
void GetPoints() {

	anim.avgEnergy = 0.f;
	int t = 0, num = 1;
	while (t < HALFN && num < HALFN) {
		dataVals.points[num][0] = dataVals.xVals[t];
		float y;

		// apply appropriate fft mode
		if (display.fftmode == 1)
			y = 0.6f + display.volumeAdjust * log10f(fft.out[t].r * fft.out[t].r + fft.out[t].i * fft.out[t].i);

		// default no fft
		else y = 64 * display.volumeAdjust * fft.inBuffer[t];

		// move each value that changed by an animation-rate adjusted amount
		float prev = dataVals.previous[num];
		if (y > prev) {
			if (display.yscaling > 0)
				dataVals.points[num][1] = prev * (1 + anim.upFactor * anim.rateAdjust);
			else dataVals.points[num][1] = y;
		}
		else {
			if (display.yscaling > 0)
				dataVals.points[num][1] = prev * (1 - anim.downfactor * anim.rateAdjust);
			else dataVals.points[num][1] = y;
		}

		// add point into avg energy sum
		anim.avgEnergy += y;

		// update previous value
		dataVals.previous[num] = dataVals.points[num][1];
		num++;
		int current = t;
		t++;

		// increase the counter until desired x increment is reached
		while (t < HALFN && dataVals.xVals[t] - dataVals.xVals[current] < display.xInc) t++;
	}

	// set first and last point to 0
	dataVals.points[0][0] = dataVals.xVals[0];
	dataVals.points[0][1] = 0;
	dataVals.points[num][0] = dataVals.xVals[HALFN - 1];
	dataVals.points[num][1] = 0;

	// store the number of points
	num++;
	if (dataVals.numPoints == 0) dataVals.numPoints = num;

	// calculate average energy
	anim.avgEnergy /= num;
	anim.avgEnergy *= anim.avgEnergy;
}

/*
* Clear out all point data and fft buffers
*/
void ClearPoints() {
	for (int i = 0; i < HALFN; i++) {
		dataVals.points[i][0] = 0;
		dataVals.points[i][1] = 0;
		dataVals.previous[i] = 0;
	}
	for (int j = 0; j < N; j++) {
		fft.inBuffer[j] = 0;
		fft.in[j].i = 0;
		fft.in[j].r = 0;
		fft.out[j].i = 0;
		fft.out[j].r = 0;
	}
	anim.peak = 0;
}

/*
* Clear all vertices, colours, texture coordinates above a given point
*/
void ClearVertices(int startAt) {
	for (int i = startAt; i < MAXVERTS; i++) {
		dataVals.vertices[i][0] = 0;
		dataVals.vertices[i][1] = 0;
		dataVals.vertices[i][2] = 0;
		dataVals.colours[i][0] = 0;
		dataVals.colours[i][1] = 0;
		dataVals.colours[i][2] = 0;
		dataVals.texCoords[i][0] = 0;
		dataVals.texCoords[i][1] = 0;
		dataVals.texCoords[i][2] = 0;
	}
}

/*
* Clear all index data above a given point
*/
void ClearIndices(int startAt) {
	for (int i = startAt; i < MAXINDS; i++)
		dataVals.indices[i] = 0;
}

/*
 * Determine and store indices for main shape
 */
void GetIndices() {

	int i = dataVals.startVertex, ii = dataVals.startIndex;
	int max = 2;//(int)floor(1.0 / bspline.uinc);

	while (i + 1 < MAXVERTS - max && ii + 2 < MAXINDS) {
		dataVals.indices[ii++] = i;
		dataVals.indices[ii++] = i + 1;
		dataVals.indices[ii++] = i + 2;
		
		i += 3;
	}
}

/*
* Animate the colour of the sky
*/
void UpdateSkyColour() {

	float r, g, b;
	int ring = display.skySphereSplits * 3;

	if (flags.colourSky) {
		float baseColour = anim.peak * anim.peak;
		float sinMod = 0.05f + 0.9f * (1 + sin(0.511f * anim.lastFrameTime));
		float cosMod = 0.1f + 0.7f * (1 + cos(0.789f * anim.lastFrameTime));
		float cosMod2 = 0.15f + 0.6f * (1 + cos(0.303f * anim.lastFrameTime));
		r = sinMod * baseColour;
		g = cosMod * baseColour;
		b = cosMod2 * baseColour;
	}
	else {
		r = 0;
		g = 0;
		b = 0;
	}

	for (int k = dataVals.startVertex - 1; k >= 0; k -= ring) {
		for (int l = 0; l < ring; l++) {
			if (k - l < 0) break;

			// set colours on first ring
			if (k - l < ring) {
				dataVals.colours[k - l][0] = 0.75f
					* dataVals.colours[k - l][0]
					+ 0.25f * r;
				dataVals.colours[k - l][1] = 0.75f
					* dataVals.colours[k - l][1]
					+ 0.25f * g;
				dataVals.colours[k - l][2] = 0.75f
					* dataVals.colours[k - l][2]
					+ 0.25f * b;
			}

			// on each frame, iterate colour out along rings,
			// taking average with previous
			else {
				dataVals.colours[k - l][0] = 0.5f
					* (dataVals.colours[k - l][0]
						+ dataVals.colours[k - ring][0]);
				dataVals.colours[k - l][1] = 0.5f
					* (dataVals.colours[k - l][1]
						+ dataVals.colours[k - ring][1]);
				dataVals.colours[k - l][2] = 0.5f
					* (dataVals.colours[k - l][2]
						+ dataVals.colours[k - ring][2]);
			}
		}
	}
}

/*
* Returns a y value for a given x, corresponding to the current shape
*/
float ShapeFunction(float x) { 

	// cos shape
	if (display.shape == 1)
		return (0.5 + anim.peak) * cosf(x * HALFPI);

	// cos shape 2
	if (display.shape == 2)
		return 2 * (0.5 + anim.peak) * (1 + cosf(x * PI));

	// cos shape 3 (multiple)
	if (display.shape == 3)
		return 2 * (0.5 + anim.peak) * (1 + cosf(3 * x * PI));

	// cos shape animated
	if (display.shape == 4)
		return (0.5 + anim.peak) * (1 + cosf(3 * x * PI))
		/ (0.5 + (1 + sinf(2 * anim.colourAdjust))
			* (1 + cosf(3 * (1 + sinf(anim.colourAdjust)) * x * PI))); // TODO animation value

	// jewels shape
	else if (display.shape == 5)
		return (0.5 + anim.peak) * cosf(x * HALFPI) * cosf(3 * x * PI) * cosf(3 * x * TWOPI);

	// ellipsoid shape
	else if (display.shape == 6)
		return 2 * (0.5 + anim.peak) * sqrtf(1 - x*x);

	// spheroid shape
	else if (display.shape == 7)
		return 5 * sqrtf(1 - x*x);

	// spheroid shape animated
	else if (display.shape == 8)
		return 3 * sqrtf(1 - x*x)
			/ (0.5f + (1 + cosf(3 
				* (1 + sinf(anim.colourAdjust)) * x * TWOPI))); // TODO animation value

	// bell shape
	else if (display.shape == 9)
		return 4 * sqrtf(1 - x*x*x);

	// default no shape
	else return 0;
}

/*
* Returns an rgb array corresponding to the inputs and the current colour scheme
*/
float *ColourFunction(float x, float v, float absy) {

	float colours[] = { 1,1,1 };
	float peakAdjust = anim.peak * anim.colourAdjust;
	if (display.colours > 4) display.colours = 0;

	// rainbow without animation
	if (display.colours == 1) {
		colours[0] = 1-(0.2 + 0.5f * (1 + sinf(1 + x + v))) + anim.peak;
		colours[2] = 0.2 + 0.5f * (1 + sinf(3.f * x + v));
		colours[1] = 0.2 + absy;
		return colours;
	}

	// shifting animated rainbow
	else if (display.colours == 2) {
		colours[0] = 0.2 + 0.5f * (1 + cosf(x + v + anim.colourAdjust));
		colours[1] = 0.2 + 0.5f * (1 + sinf(3.f * x + v + anim.colourAdjust));
		colours[2] = 0.2 + absy + sinf(anim.colourAdjust);
		return colours;
	}

	// shifting animated rainbow 2
	else if (display.colours == 3) {
		colours[0] = 1.2 - absy + sinf(2.f * anim.colourAdjust);
		colours[1] = 0.2 + 0.5f * (1 + sinf(5.f * x + v + peakAdjust));
		colours[2] = 0.2 + 0.5f * (1 + cosf(3.f * x + v + peakAdjust));
		return colours;
	}

	// shifting animated rainbow stripes
	else if (display.colours == 4) {
		colours[0] = 1.2 - absy + 0.5f * (1 + sinf(13.f * (x + v + anim.colourAdjust)));
		colours[1] = 0.2 + 0.5f * (1 + sinf(19.f * (x + v + peakAdjust)));
		colours[2] = 0.2 + 0.5f * (1 + sinf(29.f * (x + v + peakAdjust)));
		return colours;
	}

	// default white
	else return colours;
}

/*
* Returns a scaled y value based on the current y-scaling scheme
*/
float *ScaleY(float x, float y, float absy) {
	float newy = y;
	float newabsy = absy;
	float result[] = { y, absy };

	// increase peak definition and emphasize highs
	if (display.yscaling == 1) {
		if (y > 0.01f) {
			newy *= (newy + 0.25f * x) * (1.5f * newy * newy + 0.5f * (x + 0.5f));
			newabsy = 1 - abs(newy);
			newy *= newy * 2;
		}
		result[0] = newy;
		result[1] = newabsy;
		return result;
	}

	// default no scaling
	else return result;
}

/*
* Determine the vertices and colours of the main shape
*/
void GetVertices() {

	// do fft
	if (display.fftmode > 0)
		kiss_fft(fft.cfg, fft.in, fft.out);

	// get control points
	GetPoints();

	if (display.mainShape == 0) {

		int index = dataVals.startIndex;
		float peakmod = 5.33f * (1 + anim.peak / 2);

		// draw main shape
		int n = 0;
		for (double u = 0; u <= 1; u += bspline.uinc) {

			// get delta value for this u
			int d = Delta(dataVals.numPoints, u);

			// generate vertex/vertices and colour for this u
			if (dataVals.numPoints >= bspline.k) {

				// get bspline values
				dvec2 point = GetBspline(dataVals.points, u, d);
				float x = point.x;

				// scale y and add shape
				float absy = 1.f;
				float *yscaled = ScaleY(x, point.y, absy);
				float y = yscaled[0];
				absy = yscaled[1];
				if (display.shape >= 0)
					y += ShapeFunction(x);
				else {
					y += dataVals.shapePoints[n];
					n++;
				}
				y *= display.yAdjust;

				// create surface of revolution points
				float xr = peakmod * x;
				for (float v = 0.f; v < display.srevMax; v += display.srevInc) { 
					dataVals.vertices[index][0] = xr;
					dataVals.vertices[index][1] = y * cosf(v);
					dataVals.vertices[index][2] = y * sinf(v);
					float *colours = ColourFunction(x, v, absy);
					dataVals.colours[index][0] = colours[0];
					dataVals.colours[index][1] = colours[1];
					dataVals.colours[index][2] = colours[2];
					dataVals.texCoords[index][2] = 0; // 0 = no texture

					if (index < MAXVERTS - 1) index++;
					else break;
				}
			}
		}

		// update sky colour if not textured
		if (!flags.textureSky) UpdateSkyColour();

		// update peak value
		if (anim.peak > 0.01f)
			anim.peak *= (1 - 0.0001f * anim.rateAdjust);
		else anim.peak = 0.f;

		// update number of indices
		dataVals.numIndices = dataVals.startIndex + 6 * (index - dataVals.startVertex);
	}

	else {

		//int max = (int)floor(1.0 / bspline.uinc);
		int index = 0;

		for (int h = dataVals.numHistory; h > 0; h--) {

			for (int i = 0; i < HALFN; i++) {
				dataVals.history[h][i][0] = dataVals.history[h - 1][i][0];
				dataVals.history[h][i][1] = dataVals.history[h - 1][i][1];
			}
		}

		// draw main shape
		int n = 0;
		for (double u = 0; u <= 1; u += bspline.uinc) {

			// get delta value for this u
			int d = Delta(dataVals.numPoints, u);

			// generate vertex/vertices and colour for this u
			if (dataVals.numPoints >= bspline.k) {

				// get bspline values
				dvec2 point = GetBspline(dataVals.points, u, d);
				float x = point.x;

				// scale y and add shape
				float absy = 1.f;
				float *yscaled = ScaleY(x, point.y, absy);
				float y = yscaled[0];
				absy = yscaled[1];
				//if (display.shape >= 0)
				//	y += ShapeFunction(x);
				//else {
				//	y += dataVals.shapePoints[n];
				//	n++;
				//}
				y *= display.yAdjust;

				//
				dataVals.history[0][index][0] = x;
				dataVals.history[0][index][1] = y;
				index++;
			}
		}

		// update history count
		if (dataVals.numHistory < dataVals.maxHistory) dataVals.numHistory++;
		int max = index;
		index = dataVals.startVertex;

		for (int i = 0; i < dataVals.numHistory, i < HALFN; i++) {
			for (int j = 0; j < max, j < HALFN; j++) {
				float x = dataVals.history[i][j][0];
				float y = dataVals.history[i][j][1];
				dataVals.vertices[index][0] = 5 * x;
				dataVals.vertices[index][1] = y;
				dataVals.vertices[index][2] = -i;
				float *colours = ColourFunction(i, j, y);
				dataVals.colours[index][0] = colours[0];
				dataVals.colours[index][1] = colours[1];
				dataVals.colours[index][2] = colours[2];
				dataVals.texCoords[index][2] = 0; // 0 = no texture
				if (index < MAXVERTS - 1) index++;
				else break;
			}
		}

		// update sky colour if not textured
		if (!flags.textureSky) UpdateSkyColour();

		// update peak value
		if (anim.peak > 0.01f)
			anim.peak *= (1 - 0.0001f * anim.rateAdjust);
		else anim.peak = 0.f;

		// update number of indices
		dataVals.numIndices = dataVals.startIndex + 6 * (index - dataVals.startVertex);
	}
}

int GetShape(int shape, float uincmod, float xfactor, float yfactor, float xtrans, float ytrans, float ztrans, float theta, float thick, int index) {
	//int index = dataVals.startIndex;
	float peakmod = 5.33f;// *(1 + anim.peak / 2);
	float uinc = uincmod * bspline.uinc;

	// draw main shape
	int n = 0;
	for (double u = 0; u <= 1; u += uinc) {

		// get delta value for this u
		int d = Delta(dataVals.numPoints, u);

		// generate vertex/vertices and colour for this u
		if (dataVals.numPoints >= bspline.k) {

			// back
			if (shape == 0) {
				// get bspline values
				dvec2 point = GetBspline(dataVals.points, u, d);
				float x = point.x;
				float dx = uinc;
				float tempx, tempy, newx, newy;

				// scale y and add shape
				float absy = 1.f;
				float *yscaled = ScaleY(x, point.y, absy);
				float y = yscaled[0];
				absy = yscaled[1];
				//y += 4 * sqrtf(1 - x*x);
				y *= display.yAdjust;
				float dy = abs(sinf(8.f*x));

				// first point
				float xr = peakmod * x;
				tempx = xr * xfactor;
				tempy = (y + dy) * yfactor;
				if (theta != 0.f) {
					float c = cosf(theta);
					float s = sinf(theta);
					newx = tempx * c - tempy * s;
					newy = tempx * s + tempy * c;
				}
				else {
					newx = tempx;
					newy = tempy;
				}
				dataVals.vertices[index][0] = newx + xtrans;
				dataVals.vertices[index][1] = newy + ytrans;
				dataVals.vertices[index][2] = 0.5f*x + ztrans;
				float *colours = ColourFunction(x, 0, absy);
				float inner = 0.f;
				//if (abs(newx + xtrans) < 5.f && abs(newy + ytrans) < 5.f)
				//	inner = 1.f;
				inner = 1 - abs(newx + xtrans) / 50.f - abs(newy + ytrans) / 50.f;

				float cadj = 0.1f + 0.55f*anim.peak * (0.8f+sinf(u*TWOPI));
				cadj = (1 - inner) + inner*cadj;
				colours[0] *= cadj;
				colours[1] *= cadj;
				colours[2] *= cadj;
				dataVals.colours[index][0] = colours[0];
				dataVals.colours[index][1] = colours[1];
				dataVals.colours[index][2] = colours[2];
				dataVals.texCoords[index][2] = 0; // 0 = no texture
				if (index < MAXVERTS - 1) index++;
				else break;

				// second point
				xr = peakmod * (x + dx);
				tempx = xr * xfactor;
				//tempy = y * yfactor;
				if (theta != 0.f) {
					float c = cosf(theta);
					float s = sinf(theta);
					newx = tempx * c - tempy * s;
					newy = tempx * s + tempy * c;
				}
				else {
					newx = tempx;
					newy = tempy;
				}
				dataVals.vertices[index][0] = newx + xtrans;
				dataVals.vertices[index][1] = newy + ytrans;
				dataVals.vertices[index][2] = 0.5f*x + ztrans;
				//colours = ColourFunction(x + dx, 0, absy);
				dataVals.colours[index][0] = colours[0];
				dataVals.colours[index][1] = colours[1];
				dataVals.colours[index][2] = colours[2];
				dataVals.texCoords[index][2] = 0; // 0 = no texture
				if (index < MAXVERTS - 1) index++;
				else break;

				// third point
				xr = peakmod * (x + dx * 0.5f + x * dx * thick);
				tempx = xr * xfactor;
				tempy = (y + dy + dx * 10) * yfactor * thick;
				if (theta != 0.f) {
					float c = cosf(theta);
					float s = sinf(theta);
					newx = tempx * c - tempy * s;
					newy = tempx * s + tempy * c;
				}
				else {
					newx = tempx;
					newy = tempy;
				}
				dataVals.vertices[index][0] = newx + xtrans;
				dataVals.vertices[index][1] = newy + ytrans;
				dataVals.vertices[index][2] = 0.5f*x + ztrans;
				//colours = ColourFunction(x + 0.5f*dx, 0, absy);
				dataVals.colours[index][0] = colours[0];
				dataVals.colours[index][1] = colours[1];
				dataVals.colours[index][2] = colours[2];
				dataVals.texCoords[index][2] = 0; // 0 = no texture
				if (index < MAXVERTS - 1) index++;
				else break;

				// set indices
				dataVals.indices[dataVals.numIndices++] = index - 1;
				dataVals.indices[dataVals.numIndices++] = index - 2;
				dataVals.indices[dataVals.numIndices++] = index - 3;

				if (u > 0) {
					dataVals.indices[dataVals.numIndices++] = index - 4;
					dataVals.indices[dataVals.numIndices++] = index - 1;
					dataVals.indices[dataVals.numIndices++] = index - 5;

					dataVals.indices[dataVals.numIndices++] = index - 5;
					dataVals.indices[dataVals.numIndices++] = index - 1;
					dataVals.indices[dataVals.numIndices++] = index - 3;
				}
			}

			// ellipsoid
			else if (shape == 2) {
				// get bspline values
				dvec2 point = GetBspline(dataVals.points, u, d);
				float x = point.x;
				float dx = 0.5f*uinc;
				float tempx, tempy, newx, newy;

				// scale y and add shape
				float absy = 1.f;
				float *yscaled = ScaleY(x, point.y, absy);
				float y = yscaled[0];
				absy = yscaled[1];
				y += 4*sqrtf(1 - x*x);
				y *= display.yAdjust;

				// first point
				float xr = peakmod * x;
				tempx = xr * xfactor;
				tempy = y * yfactor;
				if (theta != 0.f) {
					newx = tempx * cosf(theta) - tempy * sinf(theta);
					newy = tempx * sinf(theta) + tempy * cosf(theta);
				}
				else {
					newx = tempx;
					newy = tempy;
				}
				dataVals.vertices[index][0] = newx + xtrans;
				dataVals.vertices[index][1] = newy + ytrans;
				dataVals.vertices[index][2] = abs(yfactor) + ztrans;
				float *colours = ColourFunction(x, 0, absy);
				float cadj = 0.75f + 0.75f*anim.peak;
				colours[0] *= cadj;
				colours[1] *= cadj;
				colours[2] *= cadj;
				dataVals.colours[index][0] = colours[0];
				dataVals.colours[index][1] = colours[1];
				dataVals.colours[index][2] = colours[2];
				dataVals.texCoords[index][2] = 0; // 0 = no texture
				if (index < MAXVERTS - 1) index++;
				else break;

				// second point
				xr = peakmod * (x + dx);
				tempx = xr * xfactor;
				//tempy = y * yfactor;
				if (theta != 0.f) {
					newx = tempx * cosf(theta) - tempy * sinf(theta);
					newy = tempx * sinf(theta) + tempy * cosf(theta);
				}
				else {
					newx = tempx;
					newy = tempy;
				}
				dataVals.vertices[index][0] = newx + xtrans;
				dataVals.vertices[index][1] = newy + ytrans;
				dataVals.vertices[index][2] = abs(yfactor) + ztrans;
				//colours = ColourFunction(x + dx, 0, absy);
				dataVals.colours[index][0] = colours[0];
				dataVals.colours[index][1] = colours[1];
				dataVals.colours[index][2] = colours[2];
				dataVals.texCoords[index][2] = 0; // 0 = no texture
				if (index < MAXVERTS - 1) index++;
				else break;

				// third point
				xr = peakmod * (x + dx * 0.5f + x * dx * thick);
				tempx = xr * xfactor;
				tempy = (y + dx * 10) * yfactor * thick;
				if (theta != 0.f) {
					newx = tempx * cosf(theta) - tempy * sinf(theta);
					newy = tempx * sinf(theta) + tempy * cosf(theta);
				}
				else {
					newx = tempx;
					newy = tempy;
				}
				dataVals.vertices[index][0] = newx + xtrans;
				dataVals.vertices[index][1] = newy + ytrans;
				dataVals.vertices[index][2] = abs(yfactor) + abs(200*dx*yfactor * thick) + ztrans;
				//colours = ColourFunction(x + 0.5f*dx, 0, absy);
				dataVals.colours[index][0] = colours[0];
				dataVals.colours[index][1] = colours[1];
				dataVals.colours[index][2] = colours[2];
				dataVals.texCoords[index][2] = 0; // 0 = no texture
				if (index < MAXVERTS - 1) index++;
				else break;

				// set indices
				dataVals.indices[dataVals.numIndices++] = index - 1;
				dataVals.indices[dataVals.numIndices++] = index - 2;
				dataVals.indices[dataVals.numIndices++] = index - 3;

				if (u > 0) {
					dataVals.indices[dataVals.numIndices++] = index - 4;
					dataVals.indices[dataVals.numIndices++] = index - 1;
					dataVals.indices[dataVals.numIndices++] = index - 5;

					dataVals.indices[dataVals.numIndices++] = index - 5;
					dataVals.indices[dataVals.numIndices++] = index - 1;
					dataVals.indices[dataVals.numIndices++] = index - 3;
				}
			}
		}
	}
	return index;
}

/*
* Determine the vertices and colours of the main shape
*/
void GetVerts() {

	// do fft
	if (display.fftmode > 0)
		kiss_fft(fft.cfg, fft.in, fft.out);

	// get control points
	GetPoints();
	int index = 0;
	dataVals.numIndices = 0;

	// left face
	index = GetShape(2, 1.5f, -0.3f, -0.1f, -2.f, 1.f, 0.4f, PI / 12, 1.25f, index);			// face line
	index = GetShape(2, 4.f, -0.17f, -0.1f, -4.555f, 0.95f, 0.3f, 4.5 * PI / 12, 1.8f, index);	// outer face line
	index = GetShape(2, 2.f, 0.4f, 0.1f, -2.6f, 0.4f, 0.4f, 13 * PI / 12, 1.3f, index);			// lower face line
	index = GetShape(2, 12.f, -0.02f, -0.01f, -0.4f, 0.3f, 0.8f, 7 * PI / 12, 2.f, index);		// nose
	index = GetShape(2, 2.3f, 0.35f, 0.0625f, -2.3, 2.f, 0.6f, -0.4f*(anim.peak - 0.25f), 2.f, index);	// brow
	index = GetShape(2, 2.3f, 0.4f, 0.35f, -2.7f, 2.4f, 0.f, PI / 6, 1.1f, index);				// top

	// right face
	index = GetShape(2, 1.5f, -0.3f, -0.1f, 2, 1.f, 0.4f, -PI / 12, 1.25f, index);				// face line
	index = GetShape(2, 4.f, -0.17f, -0.1f, 4.555f, 0.95f, 0.3f, -4.5 * PI / 12, 1.8f, index);	// outer face line	
	index = GetShape(2, 2.f, 0.4f, 0.1f, 2.6f, 0.4f, 0.4f, -13 * PI / 12, 1.3f, index);			// lower face line
	index = GetShape(2, 12.f, -0.02f, -0.01f, 0.4f, 0.3f, 0.8f, -7 * PI / 12, 2.f, index);		// nose
	index = GetShape(2, 2.3f, 0.35f, 0.0625f, 2.3, 2.f, 0.6f, 0.4f*(anim.peak - 0.25f), 2.f, index);	// brow
	index = GetShape(2, 2.3f, 0.4f, 0.35f, 2.7f, 2.4f, 0.f, -1.25f * PI / 6, 1.1f, index);		// top

	// left eye
	index = GetShape(2, 1.8f, 0.25f, 0.0625f, -2, 1.5f, 0.5f, 0.f, 1.5f, index);		// top lid
	index = GetShape(2, 1.8f, -0.25f, -0.0625f, -2, 1.5f, 0.5f, 0.f, 1.5f, index);		// bottom lid
	index = GetShape(2, 3.f, 0.05f, 0.0125f, -2, 1.5f, 0.6f, 0.f, 5.f, index);			// top iris
	index = GetShape(2, 3.f, -0.05f, -0.0125f, -2, 1.5f, 0.6f, 0.f, 5.f, index);		// bottom iris
	index = GetShape(2, 2.f, -0.25f, -0.0625f, -2, 1.3f, 0.4f, 0.f, 1.25f, index);		// undereye
	index = GetShape(2, 6.f, -0.085f, -0.03f, -3.65f, 1.2f, 0.25f, PI/2, 1.9f, index);	// outer accent

	// right eye
	index = GetShape(2, 1.8f, 0.25f, 0.0625f, 2, 1.5f, 0.5f, 0.f, 1.5f, index);		// top lid
	index = GetShape(2, 1.8f, -0.25f, -0.0625f, 2, 1.5f, 0.5f, 0.f, 1.5f, index);	// bottom lid
	index = GetShape(2, 3.f, 0.05f, 0.0125f, 2, 1.5f, 0.6f, 0.f, 5.f, index);		// top iris
	index = GetShape(2, 3.f, -0.05f, -0.0125f, 2, 1.5f, 0.6f, 0.f, 5.f, index);		// bottom iris
	index = GetShape(2, 2.f, -0.25f, -0.0625f, 2, 1.3f, 0.4f, 0.f, 1.25f, index);	// undereye
	index = GetShape(2, 6.f, -0.085f, -0.03f, 3.65f, 1.2f, 0.25f, -PI / 2, 1.9f, index);	// outer accent

	// third eye
	index = GetShape(2, 4.f, 0.1f, 0.0625f, 0, 2.2f, 0.3f, PI / 2, 1.25f, index);	// top lid
	index = GetShape(2, 4.f, -0.1f, -0.0625f, 0, 2.2f, 0.3f, PI / 2, 1.25f, index);	// bottom lid
	index = GetShape(2, 5.f, 0.04f, 0.012f*anim.peak, 0, 2.2f, 0.4f, -PI / 2, 5.f, index);	// top iris
	index = GetShape(2, 5.f, -0.04f, -0.012f*anim.peak, 0, 2.2f, 0.4f, -PI / 2, 5.f, index);	// bottom iris
	index = GetShape(2, 8.f, 0.03f, 0.0125f, 0, 2.425f, 0.3f, 0.f, 2.f, index);		// top accent
	index = GetShape(2, 8.f, -0.03f, -0.0125f, 0, 2.f, 0.3f, 0.f, 2.f, index);		// bottom accent
	index = GetShape(2, 4.f, -0.125f, -0.0625f, 0, 1.8f, 0.3f, 0.f, 1.25f, index);	// undereye

	// mouth
	index = GetShape(2, 2.2f, 0.9f, 0.02 + 0.07 * anim.peak, 0.f, -1.3, 0.9f, 0.f, 2.f, index);		// top lip
	index = GetShape(2, 2.3f, -0.9f, -0.01 - 0.07 * anim.peak, 0.f, -1.3, 0.9f, 0.f, 3.f, index);	// bottom lip
	index = GetShape(2, 2.f, -0.8f, -0.0625f, 0.f, -2.f, 0.8f, 0.f, 1.5f, index);					// under lip
	index = GetShape(2, 9.f, 0.03f, 0.015f, -5.5f, -1.3f, 0.6f, 4.5f * PI / 12, 2.5f, index);		// left mouth accent
	index = GetShape(2, 9.f, 0.03f, 0.015f, 5.5f, -1.3f, 0.6f, -4.5f * PI / 12, 2.5f, index);		// right mouth accent
	index = GetShape(2, 3.f, 0.275f, 0.3f, -5.f, -1.3f, 0.f, PI / 2, 1.5f, index);		// left outer mouth line
	index = GetShape(2, 3.f, 0.275f, 0.3f, 5.f, -1.3f, 0.f, -PI / 2, 1.5f, index);		// right outer mouth line
	index = GetShape(2, 2.4f, -0.92f, -0.1, 0, -2.5f, 0.4f, 0.f, 1.8f, index);			// under mouth	

	// back
	float x1 = 0.f, y1 = 0.1f;
	float x2 = 10.f, y2 = 0.1f;
	float df = PI / 8.f;
	float dr = df / 4.f;
	float adj = 0.75f + anim.peak * (1.f + 0.15f * sinf(anim.colourAdjust));
	//bool even = false;
	for (float f = 0.f; f < 2 * PI; f += df) {
		float c = cosf(f);
		float s = sinf(f);
		index = GetShape(0, 4.f, -adj*3.5f, -0.4f, x1*c-y1*s, x1*s+y1*c, -5.f, f, 1.8f, index);
		float fdr = f + dr;
		c = cosf(fdr);
		s = sinf(fdr);
		index = GetShape(0, 16.f, adj*0.8f, 0.25f, x2*c - y2*s, x2*s + y2*c, -4.f, fdr, 1.8f, index);
	}

	// update sky colour if not textured
	//if (!flags.textureSky) UpdateSkyColour();

	// update peak value
	if (anim.peak > 0.001f)
		anim.peak *= (1 - 0.0002f * anim.rateAdjust);
	else anim.peak = 0.f;

	// update number of indices
	//dataVals.numIndices = dataVals.startIndex + 3 * (index - dataVals.startVertex - 2);
}

/*
 * Use the audio input to generate the main shape
 */
void GetMainShape() {

	// Request ownership of mutex
	DWORD dwCount = 0, dwWaitResult;
	while (dwCount < 20) {
		dwWaitResult = WaitForSingleObject(ghMutex, INFINITE);  // no time-out interval
		switch (dwWaitResult) {
			// The thread got ownership of the mutex
			case WAIT_OBJECT_0:
				__try {
					GetVertices();
					dwCount++;
				}
				// Release ownership of the mutex object
				__finally {	
					if (!ReleaseMutex(ghMutex)) 
						cout << "Error: could not release mutex" << endl;
				}
				break;

			// The thread got ownership of an abandoned mutex
			// The database is in an indeterminate state
			case WAIT_ABANDONED: break;
		}
	}
}

/*
* Use the audio input to generate the main shape
*/
void GetVisual() {

	// Request ownership of mutex
	DWORD dwCount = 0, dwWaitResult;
	while (dwCount < 20) {
		dwWaitResult = WaitForSingleObject(ghMutex, INFINITE);  // no time-out interval
		switch (dwWaitResult) {
			// The thread got ownership of the mutex
		case WAIT_OBJECT_0:
			__try {
				GetVerts();
				dwCount++;
			}
			// Release ownership of the mutex object
			__finally {
				if (!ReleaseMutex(ghMutex))
					cout << "Error: could not release mutex" << endl;
			}
			break;

			// The thread got ownership of an abandoned mutex
			// The database is in an indeterminate state
		case WAIT_ABANDONED: break;
		}
	}
}

/*
* Generate an array of x values for the current x-scaling scheme
*/
void GenerateXScale() {

	// log scaling
	if (display.xscaling == 1) {
		float base = 24.f;
		float logBase = log10f(base);
		float expfirst = pow(base, 1.f);
		float explast = pow(base, 3.f);
		float xstep = float(explast - expfirst) / (HALFN - 1);

		float xi = expfirst;
		for (int t = 0; t < HALFN; t++) {
			dataVals.xVals[t] = log10f(xi) / logBase - 2.f;
			xi += xstep;
			dataVals.previous[t] = 1;
		}
	}

	// default no scaling
	else {
		for (int t = 0; t < HALFN; t++) {
			dataVals.xVals[t] = (2.f * float(t) / (HALFN - 1)) - 1.f;
			dataVals.previous[t] = 1;
		}
	}
}

/*
* Get the vertices and colours to display for drawing mode
*/
void GetDrawing() {

	// draw line in the middle (y = 0) using triangles
	dataVals.vertices[0][0] = -1.f;
	dataVals.vertices[0][1] = 0.002f;
	dataVals.vertices[1][0] = -1.f;
	dataVals.vertices[1][1] = -0.002f;
	dataVals.vertices[2][0] = 1.f;
	dataVals.vertices[2][1] = 0.002f;
	dataVals.vertices[3][0] = 1.f;
	dataVals.vertices[3][1] = -0.002f;

	for (int i = 0; i < 4; i++) {
		dataVals.vertices[i][2] = 0.f;
		dataVals.colours[i][0] = 1.f;
		dataVals.colours[i][1] = 1.f;
		dataVals.colours[i][2] = 1.f;
		dataVals.texCoords[i][2] = 0.f;
	}

	dataVals.indices[0] = 0;
	dataVals.indices[1] = 1;
	dataVals.indices[2] = 2;
	dataVals.indices[3] = 2;
	dataVals.indices[4] = 1;
	dataVals.indices[5] = 3;
	dataVals.numIndices = 6;

	// display user control points drawing
	if (display.shape == -2) {

		// show control points
		int i = 4;
		for (int j = 0; j < dataVals.shapeControl.size(); j++) {

			// determine colour
			float pr = 1.f;
			float pg = 1.f;
			float pb = 1.f;
			if (dataVals.scSelected == j) {
				pg = 0.f;
				pb = 0.f;
			}

			// determine location
			float x = dataVals.shapeControl[j].x;
			float y = dataVals.shapeControl[j].y;
			float xfix = float(glWindow.h) / glWindow.w;

			// set vertices and colours (2 triangles)
			dataVals.vertices[i][0] = x + dataVals.scRadius * xfix;
			dataVals.vertices[i][1] = y + dataVals.scRadius;
			dataVals.vertices[i][2] = 0;
			dataVals.colours[i][0] = pr;
			dataVals.colours[i][1] = pg;
			dataVals.colours[i][2] = pb;
			dataVals.texCoords[i][2] = 0;
			
			dataVals.vertices[i + 1][0] = x - dataVals.scRadius * xfix;
			dataVals.vertices[i + 1][1] = y + dataVals.scRadius;
			dataVals.vertices[i + 1][2] = 0;
			dataVals.colours[i + 1][0] = pr;
			dataVals.colours[i + 1][1] = pg;
			dataVals.colours[i + 1][2] = pb;
			dataVals.texCoords[i + 1][2] = 0;
			
			dataVals.vertices[i + 2][0] = x + dataVals.scRadius * xfix;
			dataVals.vertices[i + 2][1] = y - dataVals.scRadius;
			dataVals.vertices[i + 2][2] = 0;
			dataVals.colours[i + 2][0] = pr;
			dataVals.colours[i + 2][1] = pg;
			dataVals.colours[i + 2][2] = pb;
			dataVals.texCoords[i + 2][2] = 0;
			
			dataVals.vertices[i + 3][0] = x - dataVals.scRadius * xfix;
			dataVals.vertices[i + 3][1] = y - dataVals.scRadius;
			dataVals.vertices[i + 3][2] = 0;
			dataVals.colours[i + 3][0] = pr;
			dataVals.colours[i + 3][1] = pg;
			dataVals.colours[i + 3][2] = pb;
			dataVals.texCoords[i + 3][2] = 0;

			// set indices
			dataVals.indices[dataVals.numIndices] = i;
			dataVals.indices[dataVals.numIndices + 1] = i + 1;
			dataVals.indices[dataVals.numIndices + 2] = i + 2;
			dataVals.indices[dataVals.numIndices + 3] = i + 2;
			dataVals.indices[dataVals.numIndices + 4] = i + 1;
			dataVals.indices[dataVals.numIndices + 5] = i + 3;

			dataVals.numIndices += 6;
			i += 4;
		}

		// draw bspline
		int n = 0;
		for (double u = 0; u <= 1; u += 0.001f) {

			// get delta value for this u
			int d = Delta(dataVals.shapeControl.size(), u);

			// generate value for this u
			if (dataVals.shapeControl.size() >= bspline.k) {
				dvec2 point = GetBSplineVec(dataVals.shapeControl, u, d);

				// store shape
				dataVals.pointShape[n][0] = point.x;
				dataVals.pointShape[n][1] = point.y;
				n++;

				// draw curve
				float *colour = ColourFunction(point.y*PI, point.x*PI, 1 - abs(point.y));
				dataVals.vertices[i][0] = point.x;
				dataVals.vertices[i][1] = point.y;
				dataVals.vertices[i][2] = 0;
				dataVals.colours[i][0] = colour[0];
				dataVals.colours[i][1] = colour[1];
				dataVals.colours[i][2] = colour[2];
				dataVals.texCoords[i][2] = 0;

				dataVals.vertices[i + 1][0] = point.x + 0.01f;
				dataVals.vertices[i + 1][1] = point.y + 0.005f;
				dataVals.vertices[i + 1][2] = 0;
				dataVals.colours[i + 1][0] = colour[0];
				dataVals.colours[i + 1][1] = colour[1];
				dataVals.colours[i + 1][2] = colour[2];
				dataVals.texCoords[i + 1][2] = 0;

				if (u >= 1 - bspline.uinc) break;

				dataVals.indices[dataVals.numIndices] = i;
				dataVals.indices[dataVals.numIndices + 1] = i + 1;
				dataVals.indices[dataVals.numIndices + 2] = i + 2;
				dataVals.indices[dataVals.numIndices + 3] = i + 2;
				dataVals.indices[dataVals.numIndices + 4] = i + 1;
				dataVals.indices[dataVals.numIndices + 5] = i + 3;

				dataVals.numIndices += 6;
				i += 2;
			}
		}
		dataVals.numShapeVerts = n;
	}

	// display user freehand drawing
	else {
		int i = 4;
		for (int j = 0; j < dataVals.numShapeVerts; j++) {

            // vertices and colours
			float x = dataVals.freeShape[j][0];
			float y = dataVals.freeShape[j][1];
			float *colour = ColourFunction(y*PI, x*PI, 1-abs(y));
			dataVals.vertices[i][0] = x;
			dataVals.vertices[i][1] = y;
			dataVals.vertices[i][2] = 0.f;
			dataVals.colours[i][0] = colour[0];
			dataVals.colours[i][1] = colour[1];
			dataVals.colours[i][2] = colour[2];
			dataVals.texCoords[i][2] = 0.f;

			dataVals.vertices[i + 1][0] = x + 0.01f;
			dataVals.vertices[i + 1][1] = y + 0.005f;
			dataVals.vertices[i + 1][2] = 0.f;
			dataVals.colours[i + 1][0] = colour[0];
			dataVals.colours[i + 1][1] = colour[1];
			dataVals.colours[i + 1][2] = colour[2];
			dataVals.texCoords[i + 1][2] = 0.f;

			if (j == dataVals.numShapeVerts - 1) break;

            // indices
			dataVals.indices[dataVals.numIndices] = i;
			dataVals.indices[dataVals.numIndices + 1] = i + 1;
			dataVals.indices[dataVals.numIndices + 2] = i + 2;
			dataVals.indices[dataVals.numIndices + 3] = i + 2;
			dataVals.indices[dataVals.numIndices + 4] = i + 1;
			dataVals.indices[dataVals.numIndices + 5] = i + 3;
			dataVals.numIndices += 6;

			i += 2;
		}
	}
}

/*
* Adds a point to the freehand drawing
*/
void AddFreehandPoint(float x, float y) {

	// x difference since last point
	int dx = (int)ceil(x - glWindow.mousex);
	int xdir = 1;
	if (dx < 0) {
		dx = -dx;
		xdir = -1;
	}

	// y difference since last point
	int dy = (int)ceil(y - glWindow.mousey);
	int ydir = 1;
	if (dy < 0) {
		dy = -dy;
		ydir = -1;
	}

	// fill all x values since the last point
	for (int i = 0; i < dx; i++) {
		dataVals.freeShape[dataVals.numShapeVerts][0] = 2 * (x - xdir*i) / glWindow.w - 1;
		dataVals.freeShape[dataVals.numShapeVerts][1] = -2 * (y - ydir*i*dy/dx) / glWindow.h + 1;
		dataVals.numShapeVerts++;
	}

	// update mouse location
	glWindow.mousex = x;
	glWindow.mousey = y;
}

/*
* Adds a control point to the drawing
*/
void AddControlPoint(float x, float y) {

	dataVals.shapeControl.push_back(vec2(x, y));
	GenerateKnots(dataVals.shapeControl.size());
}

/*
* Ends draw mode and implements drawn shape in visualizer
*/
void FinalizeDrawing() {
	if (!flags.drawMode) return;
	map<float, float> m;

	// finalize from control points
	if (display.shape == -2) {

		// copy shape values to a map (to sort)
		for (size_t i = 0; i < dataVals.numShapeVerts; ++i) {
			if (abs(dataVals.pointShape[i][1]) > abs(m[dataVals.pointShape[i][0]]))
				m[dataVals.pointShape[i][0]] = dataVals.pointShape[i][1];
		}

		// copy values back
		int j = 0;
		for (const auto& x : m) {
			dataVals.pointShape[j][0] = x.first;
			dataVals.pointShape[j++][1] = x.second;
		}
		dataVals.numShapeVerts = j;

        // make b-spline
		int n = 0;
		GenerateKnots(j);

		for (double u = 0; u <= 1; u += bspline.uinc) {

			// get delta value for this u
			int d = Delta(j, u);

			// generate value for this u
			if (j >= bspline.k) {
				dvec2 point = GetBspline(dataVals.pointShape, u, d);
				dataVals.shapePoints[n] = 4 * abs(point.y);
				n++;
			}
		}
	}

	// from freehand
	else if (display.shape == -1) {
		// copy shape values to a map (to sort)
		for (size_t i = 0; i < dataVals.numShapeVerts; ++i) {
			if (abs(dataVals.freeShape[i][1]) > abs(m[dataVals.freeShape[i][0]]))
				m[dataVals.freeShape[i][0]] = dataVals.freeShape[i][1];
		}

		// copy values back
		size_t j = 0;
		for (const auto& x : m) {
			dataVals.freeShape[j][0] = x.first;
			dataVals.freeShape[j++][1] = x.second;
		}
		dataVals.numShapeVerts = j;

		// convert shape to proper dimensions using bspline
		int n = 0;
		int k = bspline.k;
		float uinc = bspline.uinc;
		SetBSpline(10, uinc);
		GenerateKnots(dataVals.numShapeVerts);

		for (double u = 0; u <= 1; u += bspline.uinc) {

			// get delta value for this u
			int d = Delta(dataVals.numShapeVerts, u);

			// generate value for this u
			if (dataVals.numShapeVerts >= bspline.k) {
				dvec2 point = GetBspline(dataVals.freeShape, u, d);
				dataVals.shapePoints[n] = 4 * abs(point.y);
				n++;
			}
		}

		// reset bspline values
		SetBSpline(k, uinc);
	}

	// end drawing mode
	flags.drawMode = false;
	ClearVertices(0);
	ClearIndices(0);
	GenerateKnots(dataVals.numPoints);
	//GenerateSkySphere();
	GetIndices();
}

/*
* Move the camera to a vector matching the given angles (polar)
*/
void MoveCamera(float p, float t) {

    // keep angles <= 2 PI
	if (camera.camT < 0) camera.camT += TWOPI;
	if (camera.camT > TWOPI) camera.camT -= TWOPI;

	// invert if upside down
	if (camera.camT >= PI) camera.inverted = -1;
	else camera.inverted = 1;

	// adjust azimuth
	camera.camP += camera.inverted * p;

	// adjust altitude
	camera.camT -= t;
}

/*
* Set the scene
*/
void SetScene(int s) {

	// clear the old data
	ClearIndices(0);
	ClearVertices(0);
	ClearPoints();

	// scene 2
	if (s == 2) {
		SetBSpline(3, bspline.uinc-0.001);
	}

	// scene 3
	else if (s == 3) {
		SetBSpline(3, bspline.uinc + 0.001);
	}

	// scene 4
	else if (s == 4) {
		SetFlags(false, false, flags.antialiasing, false, true, true, true, false, true, false, true, false);
		camera.camR = 24;
		SetBSpline(3, 0.005f);
		display.xInc = 0.1; display.yAdjust = 1.f;
		display.scene = 4; display.colours = 2;
		display.shape = 1; display.mainShape = 0;
		display.skySphereSplits = 40;
		display.srevMax = TWOPI;
		display.srevSplits = 12;
		display.srevInc = 2 * PI / display.srevSplits;
		display.tessInner = 5; display.tessOuter = 5;
		display.xscaling = 1; display.yscaling = 1;
		display.fftmode = 1;
	}

	// scene 5
	else if (s == 5) {
		SetFlags(false, false, flags.antialiasing, true, true, true, true, false, true, true, false, false);
		camera.camR = 32;
		SetBSpline(3, 0.005f);
		SetDisplayValues(display.volumeAdjust * SQRTHALFN, 0.1f, 18, 12, 3, 3, 2, 7);
		display.yAdjust = 1.f;
		display.scene = 5;
		display.mainShape = 0;
		display.srevMax = TWOPI;
		display.xscaling = 1; display.yscaling = 1;
		display.fftmode = 1;
	}

	// scene 6
	else if (s == 6) {
		SetFlags(false, true, flags.antialiasing, true, true, true, true, false, true, false, true, true);
		camera.camR = 32;
		SetBSpline(3, 0.005f);
		SetDisplayValues(display.volumeAdjust * SQRTHALFN, 0.1f, 12, 16, 3, 3, 3, 8);
		display.yAdjust = 1.f;
		display.scene = 6;
		display.mainShape = 0;
		display.srevMax = TWOPI;
		display.xscaling = 1; display.yscaling = 1;
		display.fftmode = 1;
	}

	// scene 7
	else if (s == 7) {
		SetFlags(false, false, flags.antialiasing, true, false, true, true, false, true, true, false, false);
		camera.camR = 24;
		SetBSpline(3, 0.005f);
		SetDisplayValues(display.volumeAdjust * SQRTHALFN, 0.1f, 12, 32, 5, 5, 2, 4);
		display.yAdjust = 1.f;
		display.scene = 7;
		display.mainShape = 0;
		display.srevMax = TWOPI;
		display.xscaling = 1; display.yscaling = 1;
		display.fftmode = 1;
	}

	// scene 8
	else if (s == 8) {
		SetFlags(false, false, flags.antialiasing, false, true, true, true, false, true, true, false, false);
		camera.camR = 36;
		SetBSpline(3, 0.005f);
		SetDisplayValues(display.volumeAdjust * SQRTHALFN, 0.1f, 3, 6, 2, 1, 3, 7);
		display.yAdjust = 1.f;
		display.scene = 8;
		display.mainShape = 0;
		display.srevMax = TWOPI;
		display.xscaling = 1; display.yscaling = 1;
		display.fftmode = 1;
	}

	// scene 9
	else if (s == 9) {
		SetFlags(false, false, flags.antialiasing, false, true, true, true, false, true, false, true, false);
		camera.camR = 40;
		SetBSpline(3, 0.005f);
		SetDisplayValues(display.volumeAdjust * SQRTHALFN, 0.1f, 16, 64, 6, 6, 4, 9);
		display.yAdjust = 1.f;
		display.scene = 9;
		display.mainShape = 0;
		display.srevMax = TWOPI;
		display.xscaling = 1; display.yscaling = 1;
		display.fftmode = 1;
	}

	// default to scene 1
	else {
		SetFlags(false, false, flags.antialiasing, false, false, false, false, true, false, false, false, false);
		camera.camR = 20;
		SetBSpline(2, 0.001f);
		SetDisplayValues(display.volumeAdjust * SQRTHALFN, 0.1f, 2, 6, 1, 1, 1, 0);
		display.yAdjust = 1.f;
		display.scene = 1;
		display.mainShape = 0;
		display.srevMax = display.srevInc*2;
		display.xscaling = 0; display.yscaling = 0;
		display.fftmode = 0;
	}

	// create the scene
	GenerateKnots(dataVals.numPoints);
	GenerateXScale();
	//GenerateSkySphere();
	GetIndices();
}
// --------------------------------------------------------------------------
/*
* Create buffers and fill with geometry data, returning true if successful
*/
bool InitGeometry(MyGeometry *geometry) {

    // get vertices, colours, and indices
	//if (!flags.drawMode) GetMainShape();
	//else GetDrawing();
	GetVisual();
    
    // set element count
	geometry->elementCount = dataVals.numIndices;

	// these vertex attribute indices correspond to those specified for the
	// input variables in the vertex shader
	const GLuint VERTEX_INDEX = 0;
	const GLuint COLOUR_INDEX = 1;
	const GLuint TEXTURE_INDEX = 2;

	// create an array buffer object for storing our vertices
	glGenBuffers(1, &geometry->vertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, geometry->vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(dataVals.vertices), 
		dataVals.vertices, GL_STATIC_DRAW);

	// create another one for storing our colours
	glGenBuffers(1, &geometry->colourBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, geometry->colourBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(dataVals.colours), 
		dataVals.colours, GL_STATIC_DRAW);

	// create another one for storing our texture data
	glGenBuffers(1, &geometry->textureBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, geometry->textureBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(dataVals.texCoords), dataVals.texCoords, GL_STATIC_DRAW);

	// create a vertex array object encapsulating all our vertex attributes
	glGenVertexArrays(1, &geometry->vertexArray);
	glBindVertexArray(geometry->vertexArray);

	// element array -> buffer
	glGenBuffers(1, &geometry->elementBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geometry->elementBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(dataVals.indices), 
		dataVals.indices, GL_STATIC_DRAW);

	// texture array -> buffer
	glBindBuffer(GL_ARRAY_BUFFER, geometry->textureBuffer);
	glVertexAttribPointer(TEXTURE_INDEX, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(TEXTURE_INDEX);

	// associate the position array with the vertex array object
	glBindBuffer(GL_ARRAY_BUFFER, geometry->vertexBuffer);
	glVertexAttribPointer(VERTEX_INDEX, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(VERTEX_INDEX);

	// associate the colour array with the vertex array object
	glBindBuffer(GL_ARRAY_BUFFER, geometry->colourBuffer);
	glVertexAttribPointer(COLOUR_INDEX, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(COLOUR_INDEX);

	// unbind our buffers, resetting to default state
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// check for OpenGL errors and return false if error occurred
	return !CheckGLErrors();
}

/*
* Deallocate geometry-related objects
*/
void DestroyGeometry(MyGeometry *geometry) {

	// unbind and destroy our vertex array object and associated buffers
	glBindVertexArray(0);
	glDeleteVertexArrays(1, &geometry->vertexArray);
	glDeleteBuffers(1, &geometry->vertexBuffer);
	glDeleteBuffers(1, &geometry->colourBuffer);
	glDeleteBuffers(1, &geometry->elementBuffer);
	glDeleteBuffers(1, &geometry->textureBuffer);
}

// --------------------------------------------------------------------------
/*
* Rendering function that draws our scene to the frame buffer
*/
void RenderScene(MyGeometry *geometry, MyShader *shader) {//, MyTexture *textures) {

	// toggle wireframe only
	if (flags.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// clear screen
	glClearColor(0, 0, 0, 1.0f);
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set antialiasing/blending
	if (flags.antialiasing) {
		glEnable(GL_MULTISAMPLE);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_POLYGON_SMOOTH);
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //GL_SRC_ALPHA_SATURATE, GL_ONE);
	}
	else {
		glDisable(GL_MULTISAMPLE);
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_POLYGON_SMOOTH);
		glDisable(GL_BLEND);
	}

	// bind our shader program and the vertex array object containing our
	// scene geometry, then tell OpenGL to draw our geometry
	glUseProgram(shader->program);
	glBindVertexArray(geometry->vertexArray);

	// bind textures
	for (int i = 0; i < NUMTEX; i++) {
		if (i == 0 && !flags.textureSky) continue;
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(dataVals.texture.target, dataVals.texture.textureID);
	}

	// draw patches of 3
	glPatchParameteri(GL_PATCH_VERTICES, 3);
	glDrawElements(GL_PATCHES, geometry->elementCount, GL_UNSIGNED_INT, NULL);

	// reset state to default (no shader or geometry bound)
	glBindTexture(dataVals.texture.target, 0);
	glActiveTexture(GL_TEXTURE0 + 0);
	glBindVertexArray(0);
	glUseProgram(0);

	// check for an report any OpenGL errors
	CheckGLErrors();

	// yield
	Sleep(0);
}

// --------------------------------------------------------------------------
// GLFW callback functions

/*
* Reports GLFW errors
*/
void ErrorCallback(int error, const char* description) {
	cout << "GLFW ERROR: " << error << ":" << endl;
	cout << description << endl;
}

/*
* Handles keyboard input events
*/
void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {

	// close window
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	// finalize drawing
	if (key == GLFW_KEY_ENTER && action == GLFW_PRESS)
		FinalizeDrawing();

	// toggle wireframe
	if (key == GLFW_KEY_W && action == GLFW_PRESS)
		flags.wireframe = !flags.wireframe;

	// toggle antialiasing
	if (key == GLFW_KEY_A && action == GLFW_PRESS)
		flags.antialiasing = !flags.antialiasing;

	// toggle peak rotate
	if (key == GLFW_KEY_R && action == GLFW_PRESS)
		flags.peakRotateP = !flags.peakRotateP;
	if (key == GLFW_KEY_T && action == GLFW_PRESS)
		flags.peakRotateT = !flags.peakRotateT;

	// toggle lighting
	if (key == GLFW_KEY_L && action == GLFW_PRESS)
		flags.lighting = !flags.lighting;

	// toggle simple normals
	if (key == GLFW_KEY_N && action == GLFW_PRESS)
		flags.simpleNormals = !flags.simpleNormals;

	// toggle switch outer normal
	if (key == GLFW_KEY_M && action == GLFW_PRESS) {
		flags.switchOuterNorm = !flags.switchOuterNorm;
		//GenerateSkySphere();
	}

	// toggle tess animation
	if (key == GLFW_KEY_D && action == GLFW_PRESS)
		flags.animTess = !flags.animTess;

	// toggle explode
	if (key == GLFW_KEY_E && action == GLFW_PRESS)
		flags.explode = !flags.explode;

	// toggle sky texture
	if (key == GLFW_KEY_S && action == GLFW_PRESS) {
		flags.textureSky = !flags.textureSky;
		flags.colourSky = !flags.colourSky;
	}

	// increase/decrease inner tessellation
	if (key == GLFW_KEY_HOME && action == GLFW_PRESS)
		if (display.tessInner < 11) display.tessInner++;
	if (key == GLFW_KEY_END && action == GLFW_PRESS)
		if (display.tessInner > 1) display.tessInner--;

	// increase/decrease outer tessellation
	if (key == GLFW_KEY_PAGE_UP && action == GLFW_PRESS)
		if (display.tessOuter < 11) display.tessOuter++;
	if (key == GLFW_KEY_PAGE_DOWN && action == GLFW_PRESS)
		if (display.tessOuter > 1) display.tessOuter--;

	// set scene
	if (key == GLFW_KEY_KP_0 && action == GLFW_PRESS)
		SetScene(0);
	if (key == GLFW_KEY_KP_1 && action == GLFW_PRESS)
		SetScene(1);
	if (key == GLFW_KEY_KP_2 && action == GLFW_PRESS)
		SetScene(2);
	if (key == GLFW_KEY_KP_3 && action == GLFW_PRESS)
		SetScene(3);
	if (key == GLFW_KEY_KP_4 && action == GLFW_PRESS)
		SetScene(4);
	if (key == GLFW_KEY_KP_5 && action == GLFW_PRESS)
		SetScene(5);
	if (key == GLFW_KEY_KP_6 && action == GLFW_PRESS)
		SetScene(6);
	if (key == GLFW_KEY_KP_7 && action == GLFW_PRESS)
		SetScene(7);
	if (key == GLFW_KEY_KP_8 && action == GLFW_PRESS)
		SetScene(8);
	if (key == GLFW_KEY_KP_9 && action == GLFW_PRESS)
		SetScene(9);

	// set shape
	if (key == GLFW_KEY_1 && action == GLFW_PRESS)
		display.shape = 0;
	if (key == GLFW_KEY_2 && action == GLFW_PRESS)
		display.shape = 1;
	if (key == GLFW_KEY_3 && action == GLFW_PRESS)
		display.shape = 2;
	if (key == GLFW_KEY_4 && action == GLFW_PRESS)
		display.shape = 3;
	if (key == GLFW_KEY_5 && action == GLFW_PRESS)
		display.shape = 4;
	if (key == GLFW_KEY_6 && action == GLFW_PRESS)
		display.shape = 5;
	if (key == GLFW_KEY_7 && action == GLFW_PRESS)
		display.shape = 6;
	if (key == GLFW_KEY_8 && action == GLFW_PRESS)
		display.shape = 7;
	if (key == GLFW_KEY_9 && action == GLFW_PRESS)
		display.shape = 8;
	if (key == GLFW_KEY_0 && action == GLFW_PRESS)
		display.shape = 9;
	if (key == GLFW_KEY_MINUS && action == GLFW_PRESS) {
		display.shape = -1;
		flags.drawMode = true;
	}
	if (key == GLFW_KEY_EQUAL && action == GLFW_PRESS) {
		display.shape = -2;
		flags.drawMode = true;
		GenerateKnots(dataVals.shapeControl.size());
	}

	// set fft mode
	if (key == GLFW_KEY_F && action == GLFW_PRESS)
		display.fftmode = !display.fftmode;

	// set colour mode
	if (key == GLFW_KEY_C && action == GLFW_PRESS)
		display.colours++;

	// set x scaling
	if (key == GLFW_KEY_X && action == GLFW_PRESS) {
		display.xscaling = !display.xscaling;
		GenerateXScale();
	}

	// set y scaling
	if (key == GLFW_KEY_Y && action == GLFW_PRESS)
		display.yscaling = !display.yscaling;

	// increase/decrease sky sphere splits
	if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {
		if (display.skySphereSplits < 40) display.skySphereSplits++;
		ClearVertices(0);
		ClearIndices(0);
		//GenerateSkySphere();
		GetIndices();
	}
	if (key == GLFW_KEY_LEFT && action == GLFW_PRESS) {
		if (display.skySphereSplits > 1) display.skySphereSplits--;
		ClearVertices(0);
		ClearIndices(0);
		//GenerateSkySphere();
		GetIndices();
	}

	// increase/decrease surface of revolution splits
	if (key == GLFW_KEY_UP && action == GLFW_PRESS) {
		if (display.srevSplits < 40) display.srevSplits++;
		display.srevInc = PI * 2 / display.srevSplits;
		ClearVertices(dataVals.startVertex);
		ClearIndices(dataVals.startIndex);
		GetIndices();
	}
	if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {
		if (display.srevSplits > 1) display.srevSplits--;
		display.srevInc = PI * 2 / display.srevSplits;
		ClearVertices(dataVals.startVertex);
		ClearIndices(dataVals.startIndex);
		GetIndices();
	}

	// increase/decrease volume scaling factor
	if (key == GLFW_KEY_KP_ADD && action == GLFW_PRESS)
		if (display.volumeAdjust * SQRTHALFN < 7.5f) display.volumeAdjust *= 1.1f;
	if (key == GLFW_KEY_KP_SUBTRACT && action == GLFW_PRESS) 
		if (display.volumeAdjust * SQRTHALFN > 0.5f) display.volumeAdjust *= 0.9f;
}

/*
* Handles mouse button events
*/
void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {

    // get and store mouse position
	glfwGetCursorPos(window, &glWindow.mousex, &glWindow.mousey);
    
    // convert pixels to GL coordinates
	double x = (2 * glWindow.mousex / glWindow.w) - 1;
	double y = (-2 * glWindow.mousey / glWindow.h) + 1;

	// left click
	if (button == GLFW_MOUSE_BUTTON_LEFT) {

		// when pressed start clicking
		if (action == GLFW_PRESS) {
			flags.clicking = true;

			// begin freehand drawing
			if (flags.drawMode == true && display.shape == -1) {
				dataVals.numShapeVerts = 0;
				AddFreehandPoint(glWindow.mousex, glWindow.mousey);
			}

			// control point drawing
			else if (flags.drawMode == true && display.shape == -2) {
				dataVals.scSelected = -1;

				// move a point
				for (int i = 0; i < dataVals.shapeControl.size(); i++) {
					if (abs(dataVals.shapeControl[i].x - x) <= dataVals.scRadius
						&& abs(dataVals.shapeControl[i].y - y) <= dataVals.scRadius) {
						dataVals.scSelected = i;
						flags.clicking = true;
					}
				}

				// add a point
				if (dataVals.scSelected == -1)
					AddControlPoint(x, y);
			}
		}

		// when released stop clicking
		else if (action == GLFW_RELEASE) 
			flags.clicking = false;
	}

	// right click
	else if (button == GLFW_MOUSE_BUTTON_RIGHT) {

		// if drawing control points
		if (flags.drawMode == true && display.shape == -2 && dataVals.scSelected >= 0) {

			// select a point
			for (int i = 0; i < dataVals.shapeControl.size(); i++) {
				if (abs(dataVals.shapeControl[i].x - x) <= dataVals.scRadius
					&& abs(dataVals.shapeControl[i].y - y) <= dataVals.scRadius) {
					dataVals.scSelected = i;
				}
			}

			// delete selected control point
			dataVals.shapeControl.erase(dataVals.shapeControl.begin() + dataVals.scSelected);
			dataVals.scSelected = -1;
			GenerateKnots(dataVals.shapeControl.size());
		}
	}
}

/*
* Handles cursor moving events
*/
void CursorPosCallback(GLFWwindow* window, double xpos, double ypos) {

	// if clicking, move camera based on cursor position
	if (flags.clicking) {

		// rotate camera
		if (flags.drawMode == false) {
			MoveCamera(camera.scrollSpeed * double(xpos - glWindow.mousex) / double(glWindow.w),
					camera.scrollSpeed * double(ypos - glWindow.mousey) / double(glWindow.h));

			// update current mouse position
			glWindow.mousex = xpos;
			glWindow.mousey = ypos;
		}

		// drawing control points
		else if (display.shape == -2) {

			// move a point
			if (flags.clicking && dataVals.scSelected >= 0) {
				double newx = (2 * xpos / glWindow.w) - 1;
				double newy = (-2 * ypos / glWindow.h) + 1;
				dataVals.shapeControl[dataVals.scSelected] = vec2(newx, newy);
			}
		}

		// draw freehand shape
		else 
			AddFreehandPoint(xpos, ypos);
	}
}

/*
* Handles mousewheel events
*/
void ScrollCallback(GLFWwindow* window, double xOffset, double yOffset) {

	// adjust zoom level
	camera.camR -= camera.zoomSpeed * yOffset;
	if (camera.camR < camera.minDistance) camera.camR = camera.minDistance;
	if (camera.camR > camera.maxDistance) camera.camR = camera.maxDistance;
}
// --------------------------------------------------------------------------
/*
* Initializes the GLFW window and GLAD
*/
GLFWwindow *InitGLWindow() {

	// initialize the GLFW windowing system
	if (!glfwInit()) {
		cout << "ERROR: GLFW failed to initialize" << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
	glfwSetErrorCallback(ErrorCallback);

	// set window values
	GLFWwindow *window = 0;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// set multisampling
	glfwWindowHint(GLFW_SAMPLES, glWindow.samples);

	// attempt to create a window with an OpenGL 4.1 core profile context
	window = glfwCreateWindow(glWindow.initW, glWindow.initH, "Visualizer", 0, 0);
	if (!window) {
		cout << "ERROR: failed to create GLFW window" << endl;
		cin.get();
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	// set callback functions and make our context current (active)
	glfwSetKeyCallback(window, KeyCallback);
	glfwSetMouseButtonCallback(window, MouseButtonCallback);
	glfwSetCursorPosCallback(window, CursorPosCallback);
	glfwSetScrollCallback(window, ScrollCallback);
	glfwMakeContextCurrent(window);

	// Intialize GLAD
	if (!gladLoadGL()) {
		cout << "ERROR: GLAD init failed" << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}

	// query and print out information about our OpenGL environment
	QueryGLVersion();

	return window;
}

/*
* Initializes the mutex used for audio data access
*/
void InitMutex() {
	// default security attributes, initially not owned, unnamed mutex
	ghMutex = CreateMutex(NULL, FALSE, NULL);
	if (ghMutex == NULL) {
		printf("ERROR: failed to create mutex - %d\n", GetLastError());
		cin.get();
		exit(EXIT_FAILURE);
	}
}

/*
* Initializes and starts the input audio stream
*/
void InitAudioStream(RtAudio *audio) {

	unsigned int bufferFrames = N;

	// set stream parameters
	RtAudio::StreamParameters parameters;
	parameters.deviceId = ain.device; 
	parameters.nChannels = ain.numChannels;
	parameters.firstChannel = ain.channel;

	// set stream options
	RtAudio::StreamOptions options;
	options.flags = RTAUDIO_NONINTERLEAVED || RTAUDIO_MINIMIZE_LATENCY;

	// try to open and start audio stream
	try {
		audio->openStream(NULL, &parameters, RTAUDIO_FLOAT32, ain.sampleRate,
				&bufferFrames, &AudioCallback, NULL, &options, NULL);
		audio->startStream();
	}
	catch (RtAudioError& e) {
		// Print any errors
		e.printMessage();
		cin.get();
		exit(EXIT_FAILURE);
	}
}

/*
* Terminates the input audio stream
*/
void TerminateAudioStream(RtAudio *audio) {
	try {
		// Stop the stream
		audio->stopStream();
		if (audio->isStreamOpen()) audio->closeStream();
	}
	catch (RtAudioError& e) {
		// Print any errors
		e.printMessage();
		cin.get();
	}
}

/*
* Initializes locations of uniform vars
*/
void InitUniforms(MyUniforms *uniforms, MyShader *shader) {

	glUseProgram(shader->program);
	uniforms->viewL = glGetUniformLocation(shader->program, "view");
	uniforms->projL = glGetUniformLocation(shader->program, "proj");
	uniforms->camLocL = glGetUniformLocation(shader->program, "camPoint");

	uniforms->peakL = glGetUniformLocation(shader->program, "peak");
	uniforms->simpleNormL = glGetUniformLocation(shader->program, "simpleNorm");
	uniforms->tiL = glGetUniformLocation(shader->program, "TessLevelInner");
	uniforms->toL = glGetUniformLocation(shader->program, "TessLevelOuter");
	uniforms->explodeL = glGetUniformLocation(shader->program, "explode");
	uniforms->textureSkyL = glGetUniformLocation(shader->program, "textureSky");
	uniforms->drawL = glGetUniformLocation(shader->program, "draw");

	uniforms->lightingL = glGetUniformLocation(shader->program, "showLight");
	uniforms->lightLocL = glGetUniformLocation(shader->program, "light");
	uniforms->ambientL = glGetUniformLocation(shader->program, "ambient");
	uniforms->diffRatioL = glGetUniformLocation(shader->program, "diffRatio");
	uniforms->intensityL = glGetUniformLocation(shader->program, "intensity");
	uniforms->phongL = glGetUniformLocation(shader->program, "phong");
	uniforms->specColourL = glGetUniformLocation(shader->program, "specColour");

	// set texture uniforms
	glUniform1i(glGetUniformLocation(shader->program, "tex0"), 0);
}

/*
* Sends values to uniform vars
*/
void UpdateUniforms(MyUniforms *uniforms, MyShader *shader) {

	glUseProgram(shader->program);
	glUniform1f(uniforms->peakL, anim.peak);
	glUniform1i(uniforms->simpleNormL, flags.simpleNormals);
	glUniform1f(uniforms->tiL, display.tessInner);
	glUniform1f(uniforms->toL, display.tessOuter);
	glUniform1i(uniforms->explodeL, flags.explode);
	glUniform1i(uniforms->textureSkyL, flags.textureSky);
	glUniform1i(uniforms->drawL, flags.drawMode);

	glUniformMatrix4fv(uniforms->viewL, 1, false, value_ptr(uniforms->viewMatrix));
	glUniformMatrix4fv(uniforms->projL, 1, false, value_ptr(uniforms->projMatrix));
	glUniform3fv(uniforms->camLocL, 1, uniforms->camPoint);

	glUniform1i(uniforms->lightingL, flags.lighting);
	glUniform3fv(uniforms->lightLocL, 1, light.loc);
	glUniform1f(uniforms->ambientL, light.ambient);
	glUniform1f(uniforms->diffRatioL, light.diffRatio);
	glUniform1f(uniforms->intensityL, light.intensity);
	glUniform1f(uniforms->phongL, light.phong);
	glUniform3fv(uniforms->specColourL, 1, light.specColour);
}

/*
* Generates view and projection matrixes to send to shaders
*/
void UpdateCamera(MyUniforms *uniforms) {

	// axes and translation vectors
	//vec3 xaxis = vec3(1, 0, 0);
	vec3 yaxis = vec3(0, 1, 0);
	//vec3 zaxis = vec3(0, 0, 1);

	// determine aspect ratio
	float aspectRatio = float(glWindow.w) / glWindow.h;

	// spherical camera
	float camX = camera.camR * cos(camera.camP) * sin(camera.camT);
	float camY = camera.camR * cos(camera.camT);
	float camZ = camera.camR * sin(camera.camP) * sin(camera.camT);
	vec3 cameraLoc(camX, camY, camZ);
	vec3 cameraDir = vec3(0.0, 0.0, 0.0) - cameraLoc;
	vec3 cx = cross(yaxis, cameraDir);
	vec3 cameraUp = camera.inverted * normalize(cross(cameraDir, cx));

	// field of view animation
	float fov = camera.fov + 0.05f*anim.peak * sin(anim.colourAdjust) * PI / 36; // TODO smooth peak, animation val

	// set uniform values
	uniforms->viewMatrix = lookAt(cameraLoc, cameraLoc + cameraDir, cameraUp);
	uniforms->projMatrix = perspective(fov, aspectRatio, camera.zNear, camera.zFar);
	uniforms->camPoint[0] = camX;
	uniforms->camPoint[1] = camY;
	uniforms->camPoint[2] = camZ;
}

/*
* Updates animation values each frame
*/
void UpdateAnimations() {

	// set animation speed adjustment factor
	anim.rateAdjust = (glfwGetTime() - anim.lastFrameTime) * anim.targetFps;

	// colour animation
	if (anim.colourAdjust > TWOPI) anim.colourAdjust -= TWOPI;
	anim.colourAdjust += anim.colourSpeed * anim.rateAdjust;

	// autorotation
	if (flags.peakRotateP)
		MoveCamera(0.0027f*cosf(anim.colourAdjust) * anim.peak, 0);
	if (flags.peakRotateT) 
		MoveCamera(0, 0.0027f*sinf(anim.colourAdjust) * anim.peak);

	// tessellation animation
	if (flags.animTess) {
		/* old animation
		if (anim.avgEnergy > 0.42) {
			if (dataVals.randomIndex + 1 >= HALFN) 
				dataVals.randomIndex = 0;
			display.tessInner = dataVals.randoms[dataVals.randomIndex];
			display.tessOuter = dataVals.randoms[dataVals.randomIndex + 1];
			dataVals.randomIndex += 2;
		}*/
		//int tess = 1 + int(anim.avgEnergy * 6.f);
		//if (display.tessInner > tess) display.tessInner--;
		//else if (display.tessInner < tess) display.tessInner++;
		//if (display.tessOuter > tess) display.tessOuter--;
		//else if (display.tessOuter < tess) display.tessOuter++;
	}

	// update time
	anim.lastFrameTime = glfwGetTime();
}

// ==========================================================================
/*
* Visualizer entry point
*/
int RunVisualizer() {

	// create structs
	RtAudio audio;
	MyUniforms uniforms;
	MyGeometry geometry;
	MyShader shader;

	// initialize mutex, audio stream, GL window, fft config
	InitMutex();
	InitAudioStream(&audio);
	glWindow.window = InitGLWindow();
	fft.cfg = kiss_fft_alloc(N, 0, 0, 0);
	
	// call function to load and compile shader programs
	if (!InitShaders(&shader, true)) {
		cout << "Error: Program could not initialize shaders" << endl;
		cin.get();
		return -1;
	}

	// initialize textures
	if (!InitTexture(&dataVals.texture, "texture.png", GL_TEXTURE_2D))
		cout << "Error: Program failed to intialize texture!" << endl;

	// initialize uniform vars and points/vertices/indices
	InitUniforms(&uniforms, &shader);
	GenerateXScale();
	GetPoints();
	GenerateKnots(dataVals.numPoints);
	//GenerateSkySphere();
	//InitRandoms(2, 7, 6, 10);
	GetIndices();

	// run an event-triggered main loop
	while (!glfwWindowShouldClose(glWindow.window)) {

		// get width and height
		glfwGetFramebufferSize(glWindow.window, &glWindow.w, &glWindow.h);
		glViewport(0, 0, glWindow.w, glWindow.h);

		// set up geometry
		if (!InitGeometry(&geometry))
			cout << "Error: failed to intialize geometry!" << endl;

		// update camera values and uniform vars
		UpdateCamera(&uniforms);
		UpdateUniforms(&uniforms, &shader);

		// draw the scene
		RenderScene(&geometry, &shader);

		// handle animations
		UpdateAnimations();

		// swap buffers, check for events, and clean up geometry
		glfwSwapBuffers(glWindow.window);
		glfwPollEvents();
		DestroyGeometry(&geometry);
	}

	// clean up and exit
	TerminateAudioStream(&audio);
	free(fft.cfg);
	CloseHandle(ghMutex);
	DestroyShaders(&shader);
	for (int i = 0; i < NUMTEX; i++)
		DestroyTexture(&dataVals.texture);
	glfwDestroyWindow(glWindow.window);
	glfwTerminate();
	
	return 0;
}

#ifdef VISUALIZER_CONSOLE

/*
* Main program entry point for console release
*/
int main() {

    // get device number
	cout << endl << "Enter a device number: (default 0)   ";
	int device;
	cin >> device;

    // get channel number
	cout << endl << "Enter a channel: (default 6)   ";
	int channel;
	cin >> channel;

    // get scene number
	cout << endl << "Enter a scene number: (default 1)   ";
	int scene;
	cin >> scene;

    // launch visualizer
	SetAudioInput(44100, 2, channel, device);
	SetScene(scene);
	RunVisualizer();
}

#endif