#include <GL/glew.h>

#include "GLReconWidget.h"
#include "UIconstants.h"

#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

#include <stdio.h>

//MSH: Ripped from cutil.h to remove dependency, replace
#  define CUDA_SAFE_CALL_NO_SYNC( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);                                            \


GLReconWidget::GLReconWidget(QWidget *parent) : QGLWidget(parent)
{
  cudaWidget = new cuGLReconWidget( MATRIX_SIZE_INITIAL_VALUE, MATRIX_SIZE_INITIAL_VALUE );
}

void GLReconWidget::setMatrixSize( unsigned int width, unsigned int height )
{
  cudaWidget->width = width;
  cudaWidget->height = height;
  
  cudaWidget->initializePBO();
}

void GLReconWidget::initializeGL()
{
  glewInit();
  
  if (!glewIsSupported("GL_VERSION_2_0 GL_VERSION_1_5 GL_ARB_vertex_buffer_object GL_ARB_pixel_buffer_object")) {
    fprintf(stderr, "Required OpenGL extensions missing.");
    exit(1);
  }
  
  cudaWidget->initializePBO();
}

void GLReconWidget::paintGL()
{
  cudaWidget->display();
}

void GLReconWidget::resizeGL( int w, int h )
{
  glViewport(0, 0, w, h);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 1.0, 0.0, 1.0, 0.0, 1.0); 
}

void GLReconWidget::mapPBO()
{
  cudaWidget->mapPBO();
}

void GLReconWidget::unmapPBO()
{
  cudaWidget->unmapPBO();
}

float* GLReconWidget::getDevPtr()
{
  return cudaWidget->getDevPtr();
}

// shader for displaying floating-point texture
static const char *shader_code = 
  "!!ARBfp1.0\n"
  "TEX result.color, fragment.texcoord, texture[0], 2D; \n"
  "END";

cuGLReconWidget::cuGLReconWidget( unsigned int width, unsigned int height )
{
  this->width = width;
  this->height = height;
  imageDevPtr = 0x0;
  pbo = texid = shader = 0;
}

GLuint cuGLReconWidget::compileASMShader(GLenum program_type, const char *code)
{
  GLuint program_id;
  glGenProgramsARB(1, &program_id);
  glBindProgramARB(program_type, program_id);
  glProgramStringARB(program_type, GL_PROGRAM_FORMAT_ASCII_ARB, (GLsizei) strlen(code), (GLubyte *) code);

  GLint error_pos;
  glGetIntegerv(GL_PROGRAM_ERROR_POSITION_ARB, &error_pos);
  if (error_pos != -1) {
    const GLubyte *error_string;
    error_string = glGetString(GL_PROGRAM_ERROR_STRING_ARB);
    fprintf(stderr, "Program error at position: %d\n%s\n", (int)error_pos, error_string);
    return 0;
  }
  return program_id;
}

void cuGLReconWidget::initializePBO()
{
  while( glGetError() != GL_NO_ERROR ){
    printf("\nWARNING: glError detected prior to initialisePBO");
    fflush(stdout);
  }
  
  // Create pixel buffer object (PBO) to "render Cuda memory" through a texture
  glGenBuffersARB(1, &pbo);
  glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);

  // Initialize PBO with zero image
  float *tmp = (float*) calloc( width*height, sizeof(float) );
  glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, width*height*sizeof(float), tmp, GL_STREAM_DRAW_ARB);

  glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
  CUDA_SAFE_CALL(cudaGLRegisterBufferObject(pbo));

  // Create texture for display
  glGenTextures(1, &texid);
  glBindTexture(GL_TEXTURE_2D, texid);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE32F_ARB, width, height, 0, GL_LUMINANCE, GL_FLOAT, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    
  glBindTexture(GL_TEXTURE_2D, 0);

  // Load shader program
  shader = compileASMShader(GL_FRAGMENT_PROGRAM_ARB, shader_code);

  while( glGetError() != GL_NO_ERROR ){
    printf("\nWARNING: glError detected prior to initialiseOpenGL");
    fflush(stdout);
  }

  free(tmp);
}

void cuGLReconWidget::mapPBO()
{
  if( width==0 || height == 0 ){
    printf("\nWARNING: pbo buffer size is 0! Has initialiseOpenGL() been called?\n");
  }

  imageDevPtr = 0x0;
  
  // Map the PBO used for rendering to Cuda device memory
  CUDA_SAFE_CALL(cudaGLMapBufferObject((void**)&imageDevPtr, pbo));
  
  if( !imageDevPtr ){
    printf("\nWARNING: no pbo allocated for reconstruction result!\n");
  }
  
  // Error check
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    printf("\nCuda error detected: %s\n", cudaGetErrorString(err) ); fflush(stdout);
    exit(1);
  }
}

void cuGLReconWidget::unmapPBO()
{
  // Unmap Cuda <-> PBO relation
  CUDA_SAFE_CALL(cudaGLUnmapBufferObject(pbo));

  // Error check
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    printf("\nCuda error detected: %s\n", cudaGetErrorString(err) ); fflush(stdout);
    exit(1);
  }
}

float* cuGLReconWidget::getDevPtr()
{
  return imageDevPtr;     
}

void cuGLReconWidget::display()
{
  // Clear window
  glClear(GL_COLOR_BUFFER_BIT);

  // Load texture from PBO
  glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);
  glBindTexture(GL_TEXTURE_2D, texid);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_LUMINANCE, GL_FLOAT, 0);
  glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

  // Use simple fragment program to display the floating point texture
  glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, shader);
  glEnable(GL_FRAGMENT_PROGRAM_ARB);
  glDisable(GL_DEPTH_TEST);

  // Render quad
  glBegin(GL_QUADS);
  {
    glVertex2f(0, 1); glTexCoord2f(0, 1);
    glVertex2f(0, 0); glTexCoord2f(0, 0);
    glVertex2f(1, 0); glTexCoord2f(1, 0);
    glVertex2f(1, 1); glTexCoord2f(1, 1);
  }
  glEnd();

  // Restore original state
  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_FRAGMENT_PROGRAM_ARB);
}
