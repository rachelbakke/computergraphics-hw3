﻿#include "pipelineProgram.h"

class BasicPipelineProgram : public PipelineProgram
{
public:
  int Init(const char *shaderBasePath);     // Init the program; "shaderBasePath" is the path to the folder containing the shaders.
  void SetModelViewMatrix(const float *m);  // m is column-major
  void SetProjectionMatrix(const float *m); // m is column-major
  void SetMode(int mode) { glUniform1i(h_mode, mode); }

protected:
  virtual int SetShaderVariableHandles();

  GLint h_projectionMatrix; // handle to the projectionMatrix variable in the shader
  GLint h_modelViewMatrix;  // handle to the modelViewMatrix variable in the shader
  GLint h_mode;
  // Note: we use the following naming convention: h_name is a handle to the shader variable "name".
};
