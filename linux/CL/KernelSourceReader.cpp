#include "KernelSourceReader.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

cl_program createProgram(cl_context context, cl_device_id device, string* fileName, cl_uint numFiles)
{
	cl_int errNum;
	cl_program program;
	string srcStdStr = "";

	for (cl_uint i = 0; i < numFiles; i++)
	{

		ifstream kernelFile(fileName[i], ios::in);
		if (!kernelFile.is_open())
		{
			cerr << "Failed to open file for reading: " << fileName <<
				endl;
			return NULL;
		}
		ostringstream oss;
		oss << kernelFile.rdbuf();
		srcStdStr += oss.str();
	}	

	const char *srcStr = srcStdStr.c_str();
	program = clCreateProgramWithSource(context, 1,
		(const char**)&srcStr,
		NULL, NULL);
	if (program == NULL)
	{
		cerr << "Failed to create CL program from source." << endl;
		return NULL;
	}
	//errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	errNum = clBuildProgram(program, 0, NULL, "-cl-mad-enable -cl-fast-relaxed-math", NULL, NULL);
	if (errNum != CL_SUCCESS)
	{
		// Determine the reason for the error
		char buildLog[16384];
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
		cerr << "Error in kernel: " << endl;
		cerr << buildLog;
		clReleaseProgram(program);
		return NULL;
	}
	return program;
}