
#ifndef VISUALIZER_H
#define VISUALIZER_H

// exported function prototypes
int RunVisualizer();
void SetGLWindow(int initW, int initH, int samples);
void SetAudioInput(int sampleRate, int numChannels, int channel, int device);
void SetDisplayValues(float volumeAdjust, float xInc, int srevSplits,
						int ssSplits, float ti, float to,
						int colour, int shape);
void SetBSpline(int k, float uinc);
void SetCamera(float scrollSpeed, float zoomSpeed, float fovDegrees);
void SetLight(float intensity, float ambient, float diffRatio,
				float phong, float specColour[3]);
void SetFlags(bool draw, bool animTess, bool aa, bool explode, bool lighting,
				bool peakFov, bool peakRotate, bool simpleNorms,
				bool switchNorms, bool texSky, bool cSky, bool wireframe);
void ClearPoints();
void ClearVertices(int startAt);
void ClearIndices(int startAt);

#endif