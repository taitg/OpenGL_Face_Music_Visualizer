
#ifndef GLSUPPORT_H
#define GLSUPPORT_H

// exported function prototypes
void QueryGLVersion();
bool CheckGLErrors();
string LoadSource(const string &filename);
GLuint CompileShader(GLenum shaderType, const string &source);
GLuint LinkProgram(GLuint vertexShader, GLuint fragmentShader);
GLuint LinkProgramTess(GLuint vertexShader, GLuint TCSshader, 
	GLuint TESshader, GLuint geometryShader, GLuint fragmentShader);
bool InitShaders(MyShader *shader, bool tess);
void DestroyShaders(MyShader *shader);
bool InitTexture(MyTexture* texture, const char* filename, 
	GLuint target);
void DestroyTexture(MyTexture *texture);

#endif
