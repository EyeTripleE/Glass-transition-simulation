#include "Initialize.h"
#include <ctime>
#include <iostream>
#include <random>

void initializeDefault(Parameters &params)
{
	params.timestep = 0.005;
	//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	cl_double targetDensity = 1.21;
	params.numParticlesType1 = 1024;
	params.numParticlesType2 = 0;
	params.totalParticles = params.numParticlesType1 + params.numParticlesType2;
	params.boundaries.s[0] = 0; //leftBoundary;
	params.boundaries.s[1] = sqrt(params.totalParticles*targetDensity);//100.792; //rightBoundary;
	params.boundaries.s[2] = 0;  //bottomBoundary;
	params.boundaries.s[3] = sqrt(params.totalParticles*targetDensity);  //topBoundary;
	params.origBoundaries.s[0] = params.boundaries.s[0];
	params.origBoundaries.s[1] = params.boundaries.s[1];
	params.oldBoundaries.s[0] = params.boundaries.s[0];
	params.oldBoundaries.s[1] = params.boundaries.s[1];
	params.simTime = 5.0; //Can be arbitrarily long or short
	params.currentTime = 0;
	params.binSize = 0.1;
	cl_double massPerUnitLength = 2;
	params.wallMass = massPerUnitLength*params.boundaries.s[3];
	cl_double forcePerUnitLength = 20;
	params.wallExternalForce = forcePerUnitLength*params.boundaries.s[3];
	/*<<<<<<<<<<<<<<<<<<<END SET UP PARAMETERS>>>>>>>>>>>>>>>>>>>>*/

	cl_double2 *position = new cl_double2[params.totalParticles];
	cl_double2 *velocity = new cl_double2[params.totalParticles];

	cl_double particleSpacing = 1.0f;
	cl_uint rightBound = (cl_uint)(params.boundaries.s[1] / particleSpacing);
	cl_uint upperBound = (cl_uint)(params.boundaries.s[3] / particleSpacing);

	if ((rightBound - 2)*(upperBound- 1) < params.totalParticles)
	{
		std::cerr << "Error: Simulation region too small to space particles" << std::endl;
		exit(7);
	}

	cl_uint i;
	cl_uint j;
	cl_uint **grid = new cl_uint*[rightBound];
	for (i = 0; i < rightBound; i++)
	{
		grid[i] = new cl_uint[upperBound];
		for (j = 0; j < upperBound; j++)
		{
			grid[i][j] = 0;
		}
	}

	//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	//cl_uint part = 0;
	std::default_random_engine generator((cl_uint)time(NULL));
	std::normal_distribution<cl_double> normalDistribution(0.0, 0.5);
	//std::uniform_real_distribution<cl_double> normalDistribution(-0.1f, 0.1f);

	std::uniform_int_distribution<cl_int> uniformDistributionX(1, rightBound - 1);
	std::uniform_int_distribution<cl_int> uniformDistributionY(0, upperBound);
	cl_uint x;
	cl_uint y;
	
	for (i = 0; i < params.totalParticles; i++)
	{
		x = uniformDistributionX(generator);
		y = uniformDistributionY(generator);
		while (grid[x][y] != 0)
		{
			if (x < rightBound - 1)
			{
				x++;
			}
			else
			{
				x = 1;
				y = (y < upperBound - 1) ? y + 1 : 1;
			}			
		}
		grid[x][y] = 1;
		position[i].s[0] = x*particleSpacing;
		position[i].s[1] = y* particleSpacing;
		velocity[i].s[0] = 0;// normalDistribution(generator);
		velocity[i].s[1] = 0;// normalDistribution(generator);
	}

	cl_int clStatus;
	//params.velocityBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_double2)*params.totalParticles, velocity, &clStatus);
	/*cl_double2 vSums;
	sumdouble2Buffer(params, params.velocityBuffer, vSums);
	vSums.s[0] /= params.totalParticles;
	vSums.s[1] /= params.totalParticles;*/
	//Make sure the velocity sums to zero
	cl_double vXSum = 0;
	cl_double vYSum = 0;
	for (cl_uint i = 0; i < params.totalParticles; i++)
	{
		vXSum += velocity[i].s[0];
		vYSum += velocity[i].s[1];
	}
	vXSum /= params.totalParticles;
	vYSum /= params.totalParticles;
	//Need a buffer subtraction function
	for (cl_uint i = 0; i < params.totalParticles; i++)
	{
		velocity[i].s[0] -= vXSum;
		velocity[i].s[1] -= vYSum;
	}

	//Set up some of the buffers
	params.velocityBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_double2)*params.totalParticles, velocity, &clStatus);
	params.positionBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2)*params.totalParticles, position, &clStatus);
	params.oldPositionBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_double2)*params.totalParticles, NULL, &clStatus);

	//Deleting arrays now that data has been transferred to buffers
	delete velocity;
	delete position;
	for (i = 0; i < rightBound; i++)
	{	
		delete[] grid[i];
	}
	delete[] grid;
}

void initializeFromFile(Parameters &params)
{
	//temporarily hardcoding
	params.timestep = 0.005;
	cl_double targetDensity = 1.21;
	params.numParticlesType1 = 1024 + 2048;
	params.numParticlesType2 = 0;
	params.totalParticles = params.numParticlesType1 + params.numParticlesType2;
	params.boundaries.s[2] = 0;  //bottomBoundary;
	params.boundaries.s[3] = sqrt(params.totalParticles*targetDensity);  //topBoundary;
	params.simTime = 50.0; //Can be arbitrarily long or short
	params.currentTime = 0;
	params.binSize = 0.1;

	/*<<<<<<<<<<<<<<<<<<<END SET UP PARAMETERS>>>>>>>>>>>>>>>>>>>>*/

	cl_double2 *position = new cl_double2[params.totalParticles];
	cl_double2 *oldPosition = new cl_double2[params.totalParticles];

	std::string fileName = "3000v2/simStateData_Thu Nov 06 2014 10h29m41s.txt";
	std::ifstream fin(fileName);

	for (cl_uint i = 0; i < params.totalParticles; ++i)
	{
		fin >> position[i].s[0];
		fin >> position[i].s[1];
		fin >> oldPosition[i].s[0];
		fin >> oldPosition[i].s[1];
	}
	fin >> params.boundaries.s[0];
	fin >> params.boundaries.s[1];
	fin >> params.oldBoundaries.s[0];
	fin >> params.oldBoundaries.s[1];
	fin >> params.origBoundaries.s[0];
	fin >> params.origBoundaries.s[1];
	fin.close();

	cl_double massPerUnitLength = 2;
	params.wallMass = massPerUnitLength*params.boundaries.s[3];
	cl_double forcePerUnitLength = 20;
	params.wallExternalForce = forcePerUnitLength*params.boundaries.s[3];
	
	//Adjust displacement to modify kinetic energy
	for (cl_uint i = 0; i < params.totalParticles; i++)
	{
		oldPosition[i].s[0] = position[i].s[0] - (position[i].s[0] - oldPosition[i].s[0]) / 1.14;
		oldPosition[i].s[1] = position[i].s[1] - (position[i].s[1] - oldPosition[i].s[1]) / 1.14;
		//oldPosition[i].s[0] = position[i].s[0];
		//oldPosition[i].s[1] = position[i].s[1];
	}

	cl_int clStatus;
	params.positionBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2)*params.totalParticles, position, &clStatus);
	params.oldPositionBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_double2)*params.totalParticles, oldPosition, &clStatus);
	
	//Deleting arrays now that data has been transferred to buffers					
	delete position;
	delete oldPosition;
}

void preSimSetup(Parameters &params, cl_device_id *device_list)
{
	//Set the local size
	clGetDeviceInfo(device_list[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &params.localSize, NULL);
	if (params.localSize > params.totalParticles)
	{
		params.localSize = params.totalParticles;
	}
	else if (params.totalParticles % params.localSize != 0)
	{
		std::cerr << "Error: global size is not a multiple of local size" << std::endl;
		exit(1);
	}

	cl_int clStatus;

	params.energyBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2)*params.totalParticles, NULL, &clStatus);
	params.accelerationBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_double2)*params.totalParticles, NULL, &clStatus);
	params.wallForceBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_double2)*params.totalParticles, NULL, &clStatus);
	params.innerVirialBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_double2)*params.totalParticles, NULL, &clStatus);
	params.wallsPEBuffer = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_double)*params.totalParticles, NULL, &clStatus);

	//Need both size_t and uint data types
	params.globalSize = params.totalParticles;
	params.maxTime = params.currentTime + params.simTime;
	params.wallsKE = 0;
	params.innerVirialSum = 0;
}