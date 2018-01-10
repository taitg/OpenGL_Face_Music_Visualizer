
#include <glad/glad.h>
#include <fstream>
#include <iterator>
#include <iostream>
#include <string>
#include "VisualizerDefs.h"
#include "soil/SOIL.h"

using namespace std;

// ==========================================================================
// SUPPORT FUNCTION DEFINITIONS

// --------------------------------------------------------------------------
// OpenGL utility functions

void QueryGLVersion() {

	// query opengl version and renderer information
	string version = reinterpret_cast<const char *>(glGetString(GL_VERSION));
	string glslver = reinterpret_cast<const char *>(glGetString(GL_SHADING_LANGUAGE_VERSION));
	string renderer = reinterpret_cast<const char *>(glGetString(GL_RENDERER));

	cout << "OpenGL [ " << version << " ] "
		<< "with GLSL [ " << glslver << " ] "
		<< "on renderer [ " << renderer << " ]" << endl;
}

bool CheckGLErrors() {

	bool error = false;
	for (GLenum flag = glGetError(); flag != GL_NO_ERROR; flag = glGetError()) {

		cout << "OpenGL ERROR:  ";
		switch (flag) {
		case GL_INVALID_ENUM:
			cout << "GL_INVALID_ENUM" << endl; break;
		case GL_INVALID_VALUE:
			cout << "GL_INVALID_VALUE" << endl; break;
		case GL_INVALID_OPERATION:
			cout << "GL_INVALID_OPERATION" << endl; break;
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			cout << "GL_INVALID_FRAMEBUFFER_OPERATION" << endl; break;
		case GL_OUT_OF_MEMORY:
			cout << "GL_OUT_OF_MEMORY" << endl; break;
		default:
			cout << "[unknown error code]" << endl;
		}
		error = true;
	}
	return error;
}

// --------------------------------------------------------------------------
// OpenGL shader support functions

/*
* Reads a text file with the given name into a string
*/
string LoadSource(const string &filename) {

	string source;
	ifstream input(filename.c_str());

	if (input) {
		copy(istreambuf_iterator<char>(input),
			istreambuf_iterator<char>(),
			back_inserter(source));
		input.close();
	}
	else {
		cout << "ERROR: Could not load shader source from file "
			<< filename << endl;
	}

	return source;
}

/*
* Creates and returns a shader object compiled from the given source
*/
GLuint CompileShader(GLenum shaderType, const string &source) {

	// allocate shader object name
	GLuint shaderObject = glCreateShader(shaderType);

	// try compiling the source as a shader of the given type
	const GLchar *source_ptr = source.c_str();
	glShaderSource(shaderObject, 1, &source_ptr, 0);
	glCompileShader(shaderObject);

	// retrieve compile status
	GLint status;
	glGetShaderiv(shaderObject, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE) {
		GLint length;
		glGetShaderiv(shaderObject, GL_INFO_LOG_LENGTH, &length);
		string info(length, ' ');
		glGetShaderInfoLog(shaderObject, info.length(), &length, &info[0]);
		cout << "ERROR compiling shader:" << endl << endl;
		cout << source << endl;
		cout << info << endl;
		cin.get();
	}

	return shaderObject;
}

/*
* Creates and returns a program object linked from vertex and fragment shaders
*/
GLuint LinkProgram(GLuint vertexShader, GLuint fragmentShader) {

	// allocate program object name
	GLuint programObject = glCreateProgram();

	// attach provided shader objects to this program
	if (vertexShader)   glAttachShader(programObject, vertexShader);
	if (fragmentShader) glAttachShader(programObject, fragmentShader);

	// try linking the program with given attachments
	glLinkProgram(programObject);

	// retrieve link status
	GLint status;
	glGetProgramiv(programObject, GL_LINK_STATUS, &status);
	if (status == GL_FALSE) {
		GLint length;
		glGetProgramiv(programObject, GL_INFO_LOG_LENGTH, &length);
		string info(length, ' ');
		glGetProgramInfoLog(programObject, info.length(), &length, &info[0]);
		cout << "ERROR linking shader program:" << endl;
		cout << info << endl;
		cin.get();
	}

	return programObject;
}
/*
* Creates and returns a program object linked from vertex and fragment shaders
*/
GLuint LinkProgramTess(GLuint vertexShader, GLuint TCSshader, GLuint TESshader, GLuint geometryShader, GLuint fragmentShader) {

	// allocate program object name
	GLuint programObject = glCreateProgram();

	// attach provided shader objects to this program
	if (vertexShader)   glAttachShader(programObject, vertexShader);
	if (TCSshader) glAttachShader(programObject, TCSshader);
	if (TESshader) glAttachShader(programObject, TESshader);
	if (geometryShader) glAttachShader(programObject, geometryShader);
	if (fragmentShader) glAttachShader(programObject, fragmentShader);

	// try linking the program with given attachments
	glLinkProgram(programObject);

	// retrieve link status
	GLint status;
	glGetProgramiv(programObject, GL_LINK_STATUS, &status);
	if (status == GL_FALSE) {
		GLint length;
		glGetProgramiv(programObject, GL_INFO_LOG_LENGTH, &length);
		string info(length, ' ');
		glGetProgramInfoLog(programObject, info.length(), &length, &info[0]);
		cout << "ERROR linking shader program:" << endl;
		cout << info << endl;
		cin.get();
	}

	return programObject;
}

// ----------------------------------------------------------------------
/*
* Load, compile, and link shaders, returning true if successful
*/
bool InitShaders(MyShader *shader, bool tess) {

	// load shader source from files
	string vertexSource = LoadSource("shaders/vertex.glsl");
	string TCSSource = LoadSource("shaders/tessControl.glsl");
	string TESSource = LoadSource("shaders/tessEval.glsl");
	string geometrySource = LoadSource("shaders/geometry.glsl");
	string fragmentSource = LoadSource("shaders/fragment.glsl");
	if (vertexSource.empty() || fragmentSource.empty()
		|| TCSSource.empty() || TESSource.empty()
		|| geometrySource.empty())
		return false;

	// compile shader source into shader objects
	shader->vertex = CompileShader(GL_VERTEX_SHADER, vertexSource);
	shader->fragment = CompileShader(GL_FRAGMENT_SHADER, fragmentSource);
	if (tess) {
		shader->TCS = CompileShader(GL_TESS_CONTROL_SHADER, TCSSource);
		shader->TES = CompileShader(GL_TESS_EVALUATION_SHADER, TESSource);
		shader->geometry = CompileShader(GL_GEOMETRY_SHADER, geometrySource);
		shader->program = LinkProgramTess(shader->vertex, shader->TCS, shader->TES, shader->geometry, shader->fragment);
	}
	else {
		shader->program = LinkProgram(shader->vertex, shader->fragment);
	}

	// check for OpenGL errors and return false if error occurred
	return !CheckGLErrors();
}

/*
* Deallocate shader-related objects
*/
void DestroyShaders(MyShader *shader) {

	// unbind any shader programs and destroy shader objects
	glUseProgram(0);
	glDeleteProgram(shader->program);
	glDeleteShader(shader->vertex);
	glDeleteShader(shader->fragment);
	glDeleteShader(shader->TCS);
	glDeleteShader(shader->TES);
	glDeleteShader(shader->geometry);
}

/*
* Initialize texture from file
*/
bool InitTexture(MyTexture* texture, const char* filename, GLuint target = GL_TEXTURE_2D) {
	
	int numComponents;
	unsigned char *data = SOIL_load_image(filename, &texture->width, &texture->height, &numComponents, 0);
	if (data != nullptr) {
	texture->target = target;
		glGenTextures(1, &texture->textureID);
		glBindTexture(texture->target, texture->textureID);
		GLuint format = numComponents == 3 ? GL_RGB : GL_RGBA;
		glTexImage2D(texture->target, 0, format, texture->width, texture->height, 0, format, GL_UNSIGNED_BYTE, data);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(texture->target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(texture->target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		// Clean up
		glBindTexture(texture->target, 0);
		SOIL_free_image_data(data);
		return !CheckGLErrors();
	}
	return true; //error
}

/*
* Deallocate texture-related objects
*/
void DestroyTexture(MyTexture *texture) {
	glBindTexture(texture->target, 0);
	glDeleteTextures(1, &texture->textureID);
}