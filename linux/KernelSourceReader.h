#include <string>

#ifndef KERNELSOURCEREADER_H
#define KERNELSOURCEREADER_H

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

cl_program createProgram(cl_context context, cl_device_id device, std::string* fileName, cl_uint numFiles);

#endif