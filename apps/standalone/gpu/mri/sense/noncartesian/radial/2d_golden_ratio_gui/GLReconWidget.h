#pragma once

#include "cuNDArray.h"

#if defined (WIN32)
#include <Windows.h>
#endif

#ifdef __MACH__
#import <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif //__MACH__

#include <QtOpenGL/QGLWidget>

class cuGLReconWidget
{
public:
  cuGLReconWidget( unsigned int width, unsigned int height );
  
  void initializePBO();
  void mapPBO();
  void unmapPBO();
  float* getDevPtr();
  void display();
  
  GLuint compileASMShader(GLenum program_type, const char *code);
  
  unsigned int width;
  unsigned int height;
  GLuint pbo;            // OpenGL pixel buffer object (map between Cuda and OpenGL)
  GLuint texid;          // Texture (display pbo)
  GLuint shader;         // Pixel shader for rendering of texture
  float *imageDevPtr;    // This is the "exchange buffer" between Cuda and OpenGL
};

class GLReconWidget : public QGLWidget
{
  Q_OBJECT
  
  public:
  GLReconWidget(QWidget* parent = 0);
  void setMatrixSize( unsigned int width, unsigned int height );
  void mapPBO();
  void unmapPBO();
  float* getDevPtr();

protected:
  void initializeGL();
  void paintGL();
  void resizeGL(int w, int h);
    
private:
  cuGLReconWidget *cudaWidget;
};
