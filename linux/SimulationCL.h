#ifndef SIMULATIONCL_H
#define SIMULATIONCL_H

#include <fstream>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

// Mandatory for using any feature of Xerces.
#include <xercesc/util/PlatformUtils.hpp>

// Use the Document Object Model (DOM) API
#include <xercesc/dom/DOM.hpp>

// Required for outputing a Xerces DOMDocument
// to a standard output stream (Also see: XMLFormatTarget)
#include <xercesc/framework/StdOutFormatTarget.hpp>

// Required for outputing a Xerces DOMDocument
// to the file system (Also see: XMLFormatTarget)
#include <xercesc/framework/LocalFileFormatTarget.hpp>

struct Parameters
{
	cl_mem positionBuffer;
	cl_mem velocityBuffer;
	cl_mem accelerationBuffer;
	cl_mem oldPositionBuffer;
	cl_mem energyBuffer;
	cl_mem wallForceBuffer;
	cl_mem innerVirialBuffer;
	cl_mem wallsPEBuffer;
	cl_kernel reducedoubleKernel;
	cl_kernel reducedouble2Kernel;
	cl_kernel forceKernel;
	cl_kernel diffEqKernel;
	cl_kernel periodicBoundaryKernel;
	//cl_kernel pairCorrelationKernel;
	cl_command_queue command_queue;
	cl_double timestep;
	cl_double4 boundaries;
	cl_uint totalParticles;
	cl_uint numParticlesType1;
	cl_uint numParticlesType2;
	size_t localSize;
	size_t globalSize;
	cl_double binSize;
	cl_double currentTime;
	cl_double simTime;
	cl_double maxTime;
	cl_double previousStepTime;
	cl_double wallMass;
	cl_double wallExternalForce;
	cl_double wallsKE;
	//cl_double wallsPE;
	cl_double innerVirialSum;
	cl_double2 oldBoundaries;
	cl_double2 wallAcceleration;
	cl_double2 origBoundaries;
	cl_double2 wallInternalForceSum;
	cl_context context;
};

void performEulerOperation(Parameters &params);

void performVerletOperation(Parameters &params);

void calculateAcceleration(Parameters &params);

void outputPosition(Parameters &params, FILE* positionFile, const cl_double currentTime);

void sumdouble2Buffer(Parameters &params, cl_mem &buffer, cl_double2& result);

void sumdoubleBuffer(Parameters &params, cl_mem &buffer, cl_double& result);

void wallVerlet(Parameters &params);

void wallEuler(Parameters &params);

void outputHistogram(Parameters &params, const char* time_string);

void outputData(Parameters &params, const char* time_string);

void outputSettings();

void DoOutput2File(xercesc::DOMDocument* pDOMDocument, const wchar_t * FullFilePath);

size_t get_nextpowerof2(size_t n);

#endif