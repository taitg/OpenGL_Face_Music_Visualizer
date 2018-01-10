
=======
Program
=======

Visualizer is a program which takes real-time input from ASIO audio drivers
and produces a visual representation of the audio.


===========
Compilation
===========

The program was developed and compiled using Visual Studio 2015. Project
files have been provided, but no makefile. All non-standard dependencies
have been included.


=====
Files
=====

Visualizer.cpp/.h		Contains main functions and methods for taking audio
						input and drawing the visualizer.

VisualizerDefs.h		Contains definitions of data structures

GLsupport.cpp/.h		Contains functions and methods for OpenGL functions
						such as loading and compiling shaders

visualizer101.cpp/.h	Main program function for generating GUI instance

visualizer101Dlg.cpp/.h	Contains logic and functions for GUI controls

The rest of the included files are libraries required for audio input,
GUI, Fourier transform, and texture loading

