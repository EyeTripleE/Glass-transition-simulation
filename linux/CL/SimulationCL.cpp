//NONSTANDARD C CODE IN BASICALLY EVERY FUNCTION WITH _s IN IT
#include "SimulationCL.h"
#include <ctime>
#include "KernelSourceReader.h"
#include "Initialize.h"
#include "DetermineVector.h"
#include "direct.h"
#include <sstream>

#pragma warning( disable : 4996) 
int main(int argc, char* argv[])
{
	/*CONSTANTS FOR REFERENCE THESE ARE HARDCODED*/
	//constants 0.5(sigma1 + sigma2)
	//double sigma1to1 = 1; //in units of sigma
	//double sigma2to2 = 1.4; //in units of sigma
	//double sigma1to2 = 1.2; // in units of sigma
	//double epsilon = 1;// in units of epsilon
	//double massParticle1 = 1; //in units of massParticle1
	//double massParticle2 = 2; //in units of massParticle1

	clock_t tstart = clock();	

	/*<< << << << << << << << << <BEGIN OPENCL STUFF >> >> >> >> >> >> >> >> >*/

	// Get platform and device information
	cl_platform_id* platforms = NULL;
	cl_uint num_platforms;

	//Set up the Platform
	cl_int clStatus = clGetPlatformIDs(0, NULL, &num_platforms);
	platforms = (cl_platform_id *)malloc(sizeof(cl_platform_id)*num_platforms);
	clStatus = clGetPlatformIDs(num_platforms, platforms, NULL);

	//Get the devices list and choose the device you want to run on
	cl_device_id *device_list = NULL;
	cl_uint num_devices;
	clStatus = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
	device_list = (cl_device_id *)malloc(sizeof(cl_device_id)*num_devices);
	clStatus = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, num_devices, device_list, NULL);
	//printf("%i", num_platforms);
	//return 0;

	// Create one OpenCL context for each device in the platform
	cl_context context;
	context = clCreateContext(NULL, num_devices, device_list, NULL, NULL, &clStatus);

	//Create program from source
	std::string fileNames[2];
	fileNames[0] = "DetermineVector.h";
	fileNames[1] = "kernelCode.cl";
	cl_program program = createProgram(context, device_list[0], fileNames, 2);

	//size_t item;
	//clGetDeviceInfo(device_list[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &item, NULL);
	//cout << item << endl;

	/*<<<<<<<<<<<<<<<<<<<END OPENCL STUFF>>>>>>>>>>>>>>>>>>>>>>>>>*/

	//Variables that are the same regardless of load method
	tm current_time;
	time_t cur_time = time(NULL);
	localtime_s(&current_time, &cur_time);
	char time_string[64];
	strftime(time_string, sizeof(time_string), "%a %b %d %Y %Hh%Mm%Ss", &current_time);
	std::string folderName = "PistonSim_";
	folderName.append(time_string);
	_mkdir(folderName.c_str());
	FILE* positionFile;
	FILE* vsTimeFile;
	FILE* timeAvgFile;
	
	std::string fileName = folderName;
	fileName.append("/position_");
	fileName.append(time_string);
	fileName.append(".txt");
	fopen_s(&positionFile, fileName.c_str(), "w");
	
	fileName = folderName;
	fileName.append("/vsTimePlots_");
	fileName.append(time_string);
	fileName.append(".txt");
	fopen_s(&vsTimeFile, fileName.c_str(), "w");

	fileName = folderName;
	fileName.append("/timeAvgResults_");
	fileName.append(time_string);
	fileName.append(".txt");
	fopen_s(&timeAvgFile, fileName.c_str(), "w");
	
	cl_double2 energies = { 0, 0 };
	//cl_double virialIntegral = 0;
	cl_double totalEnergy;
	cl_double totalKE;
	cl_double totalPE;
	cl_double virialAvg;
	cl_double boxWidth;
	cl_double boxHeight;
	cl_double area;
	cl_double virialPressure;
	cl_double wallInternPE;
	cl_double wallExternPE;
	cl_double totalWallsPE;
	cl_double avgKE;
	cl_double state;
	cl_double timeAvgTotalEnergy = 0;
	cl_double timeAvgTemp = 0;
	cl_double timeAvgPressure = 0;
	cl_double timeAvgMolarVol = 0;
	cl_double avgTurnAroundDist;
	cl_double tradPressure = 0;
	cl_double timeAvgTradPressure = 0;
	cl_double timeAvgPotentialEnergy = 0;
	cl_double timeAvgPartPotEnergy = 0;
	cl_double timeAvgBoxWidth = 0;
	cl_double timeAvgTotalPartEnergyNoWall = 0;
	cl_double timeAvgTotalPartEnergyWithWall = 0;
	cl_double timeAvgTotalPartEnergyHalfWall = 0;
	cl_double totalParticleEnergy;

	/*<<<<<<<<<<<<<<<<<<<SET UP PARAMETERS>>>>>>>>>>>>>>>>>>>>>>>>*/

	//Set up parameter structure
	Parameters params;
	params.context = context;
	params.command_queue = clCreateCommandQueue(context, device_list[0], 0, &clStatus);
	params.forceKernel = clCreateKernel(program, "calc_acceleration", &clStatus);
	params.periodicBoundaryKernel = clCreateKernel(program, "apply_periodic_boundary", &clStatus);
	params.reducedouble2Kernel = clCreateKernel(program, "reducedouble2", &clStatus);
	params.reducedoubleKernel = clCreateKernel(program, "reducedouble", &clStatus);
	//params.pairCorrelationKernel = clCreateKernel(program, "pairCorrelationKernel", &clStatus);
	char titleString[] = "CurrentTime TotalParticleEnergy TotalEnergy ParticlesPE ParticlesKE TotalPE TotalKE Walls-ParticlePE Walls-ExternPE WallsKE InstVirialPressure AvgVirialPressure BoxWidth AvgKE IdealState\n";
	fwrite(titleString, sizeof(char), sizeof(titleString) - 1, vsTimeFile);
	//Begin Simulation
	if (argc == 1)
	{
		params.diffEqKernel = clCreateKernel(program, "euler_kernel", &clStatus);
		initializeDefault(params);
		preSimSetup(params, device_list);

		outputPosition(params, positionFile, params.currentTime);
		performEulerOperation(params);
		params.currentTime += params.timestep;
		outputPosition(params, positionFile, params.currentTime);
		//sum the potential energy
		//sumBuffer(reduceKernel, params.potentialEnergyBuffer,params.totalParticles, params.localSize, potentialEnergy);
		//output the energy
		//previousStepTime = params.currentTime;
		params.currentTime += params.timestep;

		//Change the differential equation algorithm
		clStatus = clReleaseKernel(params.diffEqKernel);
		clStatus = clReleaseMemObject(params.velocityBuffer);
	}
	else if (argc == 2)
	{
		initializeFromFile(params);
		preSimSetup(params, device_list);
	}
	else
	{
		//cout << "Too many arguments, terminating" << endl;
		printf("Too many arguments, terminating\n");
		exit(3);
	}

	//Sections become unified here
	params.diffEqKernel = clCreateKernel(program, "verlet_kernel", &clStatus);
	boxHeight = params.boundaries.s[3] - params.boundaries.s[2];
	cl_uint counter = 0;	
	for (params.currentTime; params.currentTime < params.maxTime; params.currentTime += params.timestep)
	{		
		performVerletOperation(params);
		sumdouble2Buffer(params, params.energyBuffer, energies);

		avgKE = energies.s[1] / params.totalParticles;
		avgTurnAroundDist = 0;//2 * pow((3.092505268 /** pow(1, 12)*/) / (0.5*avgKE), 0.09090909090909090909090909090909);
		boxWidth = params.boundaries.s[1] - params.boundaries.s[0] - avgTurnAroundDist;
		area = (boxWidth)*(boxHeight);
		//may have to return force on piston
		//virialPressure = (energies.s[1] + params.innerVirialSum) / area;
		virialPressure = (energies.s[1] + 6*energies.s[0]) / area;
		tradPressure = (params.wallInternalForceSum.s[0] + params.wallInternalForceSum.s[1]) / (2 * boxHeight);
		timeAvgPressure += virialPressure*params.timestep;
		timeAvgMolarVol += area*params.timestep;
		timeAvgTemp += avgKE*params.timestep;
		timeAvgTradPressure += tradPressure*params.timestep;
		timeAvgBoxWidth += boxWidth * params.timestep;
		totalParticleEnergy = (energies.s[0] + energies.s[1]) / params.totalParticles;

		sumdoubleBuffer(params, params.wallsPEBuffer, wallInternPE);		
		wallExternPE = params.wallExternalForce*((params.boundaries.s[0] - params.origBoundaries.s[0]) + (params.origBoundaries.s[1] - params.boundaries.s[1]));
		totalWallsPE = wallInternPE - wallExternPE;
		totalPE = energies.s[0] + totalWallsPE;
		totalKE = energies.s[1] + params.wallsKE;
		totalEnergy = totalPE + totalKE;
		timeAvgTotalEnergy += totalEnergy*params.timestep;
		timeAvgPotentialEnergy += totalPE*params.timestep;
		timeAvgPartPotEnergy += energies.s[0] * params.timestep;
		timeAvgTotalPartEnergyNoWall += (energies.s[0] + energies.s[1])*params.timestep;
		timeAvgTotalPartEnergyWithWall += (totalPE + energies.s[1])*params.timestep;
		timeAvgTotalPartEnergyHalfWall += (energies.s[0] + totalWallsPE / 2 + energies.s[1])*params.timestep;

		//Prevent from outputting each step
		/*if (counter < 10)
		{
			counter++;
		}
		else
		{*/
			virialAvg = timeAvgPressure / params.currentTime;
			//Using molar volume for state.  What is the molar volume for a mixture of particles
			state = (virialPressure*area) / (avgKE * params.totalParticles);
			fprintf_s(vsTimeFile, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", params.currentTime, totalParticleEnergy, totalEnergy, energies.s[0], energies.s[1], totalPE, totalKE, wallInternPE, wallExternPE, params.wallsKE, virialPressure, virialAvg, boxWidth, avgKE, state);

			outputPosition(params, positionFile, params.currentTime);

			//counter = 0;
		//}
	}

	//UNCOMMENT THE NEXT LINE FOR THE HISTOGRAM
	//outputHistogram(params, folderName.c_str());
	outputData(params, folderName.c_str());

	clStatus = clReleaseMemObject(params.positionBuffer);
	clStatus = clReleaseMemObject(params.oldPositionBuffer);
	clStatus = clReleaseMemObject(params.accelerationBuffer);
	clStatus = clReleaseMemObject(params.energyBuffer);
	clStatus = clReleaseMemObject(params.wallForceBuffer);
	
	clStatus = clReleaseKernel(params.forceKernel);
	clStatus = clReleaseKernel(params.diffEqKernel);
	clStatus = clReleaseKernel(params.periodicBoundaryKernel);
	clStatus = clReleaseKernel(params.reducedoubleKernel);
	clStatus = clReleaseKernel(params.reducedouble2Kernel);

	clStatus = clReleaseProgram(program);
	clStatus = clReleaseCommandQueue(params.command_queue);
	clStatus = clReleaseContext(context);

	free(platforms);
	free(device_list);

	fclose(positionFile);
	fclose(vsTimeFile);
	printf("Time Avg Total Molar System Energy: %g\n", timeAvgTotalEnergy/(params.currentTime * params.totalParticles));
	printf("Time Avg Total Molar System Pot. Energy %g\n", timeAvgPotentialEnergy / (params.currentTime*params.totalParticles));
	printf("Time Avg Particle Molar Pot. Energy %g\n", timeAvgPartPotEnergy / (params.currentTime*params.totalParticles));
	printf("Time Avg Pressure: %g\n", timeAvgPressure/params.currentTime);
	printf("Time Avg Molar Volume: %g\n", timeAvgMolarVol/(params.currentTime * params.totalParticles));
	printf("Time Avg Temp: %g\n", timeAvgTemp/params.currentTime);
	printf("Time Avg Box Width: %g\n", timeAvgBoxWidth / params.currentTime);
	printf("Time Avg Total Molar Particle Energy No Walls PE %g\n", timeAvgTotalPartEnergyNoWall / (params.currentTime*params.totalParticles));
	printf("Time Avg Total Molar Particle Energy Half Walls PE %g\n", timeAvgTotalPartEnergyHalfWall / (params.currentTime*params.totalParticles));
	printf("Time Avg Total Molar Particle Energy With Walls PE %g\n", timeAvgTotalPartEnergyWithWall / (params.currentTime*params.totalParticles));
	printf("Total Particles %d\n", params.totalParticles);
	printf("\nRun Time: %g\n", (cl_double)(clock() - tstart) / CLOCKS_PER_SEC);

	fprintf_s(timeAvgFile,"Time Avg Total Molar Energy: %g\n", timeAvgTotalEnergy / (params.currentTime * params.totalParticles));
	fprintf_s(timeAvgFile, "Time Avg Total Molar Pot. Energy %g\n", timeAvgPotentialEnergy / (params.currentTime*params.totalParticles));
	fprintf_s(timeAvgFile, "Time Avg Particle Molar Pot. Energy %g\n", timeAvgPartPotEnergy / (params.currentTime*params.totalParticles));
	fprintf_s(timeAvgFile, "Time Avg Pressure: %g\n", timeAvgPressure / params.currentTime);
	fprintf_s(timeAvgFile, "Time Avg Molar Volume: %g\n", timeAvgMolarVol / (params.currentTime * params.totalParticles));
	fprintf_s(timeAvgFile, "Time Avg Temp: %g\n", timeAvgTemp / params.currentTime);
	fprintf_s(timeAvgFile, "Time Avg Box Width: %g\n", timeAvgBoxWidth / params.currentTime);
	fprintf_s(timeAvgFile, "Time Avg Total Molar Particle Energy No Walls PE %g\n", timeAvgTotalPartEnergyNoWall / (params.currentTime*params.totalParticles));
	fprintf_s(timeAvgFile, "Time Avg Total Molar Particle Energy Half Walls PE %g\n", timeAvgTotalPartEnergyHalfWall / (params.currentTime*params.totalParticles));
	fprintf_s(timeAvgFile, "Time Avg Total Molar Particle Energy With Walls PE %g\n", timeAvgTotalPartEnergyWithWall / (params.currentTime*params.totalParticles));
	fprintf_s(timeAvgFile, "Total Particles %d\n", params.totalParticles);
	fprintf_s(timeAvgFile, "\nRun Time: %g\n", (cl_double)(clock() - tstart) / CLOCKS_PER_SEC);

	fclose(timeAvgFile);
	return 0;
}
/*END OF MAIN FUNCTION*/

void performEulerOperation(Parameters &params)
{
	//Calculate Acceleration

	calculateAcceleration(params);
	
	//Do Euler Stuff
	clSetKernelArg(params.diffEqKernel, 0, sizeof(cl_double), (void*)&params.timestep);
	clSetKernelArg(params.diffEqKernel, 1, sizeof(cl_mem), (void*)&params.positionBuffer);
	clSetKernelArg(params.diffEqKernel, 2, sizeof(cl_mem), (void*)&params.velocityBuffer);
	clSetKernelArg(params.diffEqKernel, 3, sizeof(cl_mem), (void*)&params.accelerationBuffer);
	clSetKernelArg(params.diffEqKernel, 4, sizeof(cl_mem), (void*)&params.oldPositionBuffer);

	//Launch kernel
	clEnqueueNDRangeKernel(params.command_queue, params.diffEqKernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL);

	clFlush(params.command_queue);
	clFinish(params.command_queue);	

	//Move the walls
	wallEuler(params);

	//Apply Periodic Boundary
	cl_double intervals = params.boundaries.s[1] - params.boundaries.s[0];
	clSetKernelArg(params.periodicBoundaryKernel, 0, sizeof(cl_mem), (void*)&params.positionBuffer);
	clSetKernelArg(params.periodicBoundaryKernel, 1, sizeof(cl_mem), (void*)&params.oldPositionBuffer);
	clSetKernelArg(params.periodicBoundaryKernel, 2, sizeof(cl_double), (void*)&intervals);
	clSetKernelArg(params.periodicBoundaryKernel, 3, sizeof(cl_double4), (void*)&params.boundaries);

	clEnqueueNDRangeKernel(params.command_queue, params.periodicBoundaryKernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL);
	
	clFlush(params.command_queue);
	clFinish(params.command_queue);
}

void performVerletOperation(Parameters &params)
{
	//Calculate Acceleration
	calculateAcceleration(params);

	//Do Verlet Stuff
	clSetKernelArg(params.diffEqKernel, 0, sizeof(cl_double), (void*)&params.timestep);
	clSetKernelArg(params.diffEqKernel, 1, sizeof(cl_mem), (void*)&params.positionBuffer);
	clSetKernelArg(params.diffEqKernel, 2, sizeof(cl_mem), (void*)&params.accelerationBuffer);
	clSetKernelArg(params.diffEqKernel, 3, sizeof(cl_mem), (void*)&params.oldPositionBuffer);
	clSetKernelArg(params.diffEqKernel, 4, sizeof(cl_mem), (void*)&params.energyBuffer);
	clSetKernelArg(params.diffEqKernel, 5, sizeof(cl_uint), (void*)&params.numParticlesType1);

	//Launch kernel
	clEnqueueNDRangeKernel(params.command_queue, params.diffEqKernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL);

	clFlush(params.command_queue);
	clFinish(params.command_queue);

	//Move walls
	wallVerlet(params);

	//Apply Periodic Boundary
	cl_double intervals = params.boundaries.s[1] - params.boundaries.s[0];
	clSetKernelArg(params.periodicBoundaryKernel, 0, sizeof(cl_mem), (void*)&params.positionBuffer);
	clSetKernelArg(params.periodicBoundaryKernel, 1, sizeof(cl_mem), (void*)&params.oldPositionBuffer);
	clSetKernelArg(params.periodicBoundaryKernel, 2, sizeof(cl_double), (void*)&intervals);
	clSetKernelArg(params.periodicBoundaryKernel, 3, sizeof(cl_double4), (void*)&params.boundaries);

	clEnqueueNDRangeKernel(params.command_queue, params.periodicBoundaryKernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL);
	
	clFlush(params.command_queue);
	clFinish(params.command_queue);
}

void calculateAcceleration(Parameters &params)
{
	//Set kernel arguments
	clSetKernelArg(params.forceKernel, 0, sizeof(cl_mem), (void*) &params.positionBuffer);
	clSetKernelArg(params.forceKernel, 1, sizeof(cl_mem), (void*) &params.accelerationBuffer);
	clSetKernelArg(params.forceKernel, 2, sizeof(cl_mem), (void*) &params.wallForceBuffer);
	clSetKernelArg(params.forceKernel, 3, sizeof(cl_double4), (void*) &params.boundaries);
	clSetKernelArg(params.forceKernel, 4, sizeof(cl_uint), (void*) &params.totalParticles);
	clSetKernelArg(params.forceKernel, 5, sizeof(cl_uint), (void*)&params.numParticlesType1);
	clSetKernelArg(params.forceKernel, 6, sizeof(cl_double2)*params.localSize, NULL);
	clSetKernelArg(params.forceKernel, 7, sizeof(cl_mem), (void*) &params.energyBuffer);
	clSetKernelArg(params.forceKernel, 8, sizeof(cl_mem), (void*)&params.innerVirialBuffer);
	clSetKernelArg(params.forceKernel, 9, sizeof(cl_mem), (void*)&params.wallsPEBuffer);
	
	//Launch kernel
	clEnqueueNDRangeKernel(params.command_queue, params.forceKernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL);
	
	clFlush(params.command_queue);
	clFinish(params.command_queue);

	//cl_double2 wallForceSum;
	sumdouble2Buffer(params, params.wallForceBuffer,params.wallInternalForceSum);
	params.wallAcceleration.s[0] = (params.wallInternalForceSum.s[0] + params.wallExternalForce) / params.wallMass;
	params.wallAcceleration.s[1] = (params.wallInternalForceSum.s[1] - params.wallExternalForce) / params.wallMass;

	//wallForceSum.s[0] += params.wallExternalForce;
	//wallForceSum.s[1] -= params.wallExternalForce;
	//params.wallAcceleration.s[0] = wallForceSum.s[0] / params.wallMass;
	//params.wallAcceleration.s[1] = wallForceSum.s[1] / params.wallMass;

	cl_double tempVirialSum;
	sumdoubleBuffer(params, params.innerVirialBuffer, tempVirialSum);
	params.innerVirialSum = tempVirialSum / 4;
}

void sumdouble2Buffer(Parameters &params, cl_mem &buffer, cl_double2& result)
{	
	cl_mem resultBuffer = clCreateBuffer(params.context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2), NULL, NULL);
	
	//Set kernel arguments
	clSetKernelArg(params.reducedouble2Kernel, 0, sizeof(cl_mem), (void*)&buffer);
	clSetKernelArg(params.reducedouble2Kernel, 1, sizeof(cl_mem), (void*)&resultBuffer);	
	clSetKernelArg(params.reducedouble2Kernel, 2, sizeof(cl_double2)*get_nextpowerof2(params.localSize), NULL);
	
	//Launch kernel
	clEnqueueNDRangeKernel(params.command_queue, params.reducedouble2Kernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL);

	clFlush(params.command_queue);
	clFinish(params.command_queue);

	cl_double2* pointer = (cl_double2*)clEnqueueMapBuffer(params.command_queue, resultBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(cl_double2), 0, NULL, NULL, NULL);

	result = *pointer;

	clEnqueueUnmapMemObject(params.command_queue, resultBuffer, pointer, NULL, NULL, NULL);

	clReleaseMemObject(resultBuffer);
}

void sumdoubleBuffer(Parameters &params, cl_mem &buffer, cl_double& result)
{
	cl_mem resultBuffer = clCreateBuffer(params.context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double), NULL, NULL);

	//Set kernel arguments
	clSetKernelArg(params.reducedoubleKernel, 0, sizeof(cl_mem), (void*)&buffer);
	clSetKernelArg(params.reducedoubleKernel, 1, sizeof(cl_mem), (void*)&resultBuffer);
	clSetKernelArg(params.reducedoubleKernel, 2, sizeof(cl_double2)*get_nextpowerof2(params.localSize), NULL);

	//Launch kernel
	clEnqueueNDRangeKernel(params.command_queue, params.reducedoubleKernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL);

	clFlush(params.command_queue);
	clFinish(params.command_queue);

	cl_double* pointer = (cl_double*)clEnqueueMapBuffer(params.command_queue, resultBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(cl_double), 0, NULL, NULL, NULL);

	result = *pointer;

	clEnqueueUnmapMemObject(params.command_queue, resultBuffer, pointer, NULL, NULL, NULL);

	clReleaseMemObject(resultBuffer);
}

void wallEuler(Parameters &params)
{
	params.oldBoundaries.s[0] = params.boundaries.s[0];
	params.oldBoundaries.s[1] = params.boundaries.s[1];
	params.boundaries.s[0] += 0.5f*params.timestep*params.timestep*params.wallAcceleration.s[0];
	params.boundaries.s[1] += 0.5f*params.timestep*params.timestep*params.wallAcceleration.s[1];
}

void wallVerlet(Parameters &params)
{
	cl_double lDisp = params.boundaries.s[0] - params.oldBoundaries.s[0];
	cl_double rDisp = params.boundaries.s[1] - params.oldBoundaries.s[1];
	cl_double timestepSq = params.timestep*params.timestep;
	lDisp += timestepSq*params.wallAcceleration.s[0];
	rDisp += timestepSq*params.wallAcceleration.s[1];
	cl_double oldOldLBound = params.oldBoundaries.s[0];
	cl_double oldOldRBound = params.oldBoundaries.s[1];	
	params.oldBoundaries.s[0] = params.boundaries.s[0];
	params.oldBoundaries.s[1] = params.boundaries.s[1];
	params.boundaries.s[0] += lDisp;
	params.boundaries.s[1] += rDisp;

	cl_double lVel = (params.boundaries.s[0] - oldOldLBound) / (2 * params.timestep);
	cl_double rVel = (params.boundaries.s[1] - oldOldRBound) / (2 * params.timestep);

	params.wallsKE = 0.5f*params.wallMass*(lVel*lVel + rVel*rVel);
}

void outputPosition(Parameters &params, FILE* positionFile, const cl_double currentTime)
{
	cl_double2* pointer = (cl_double2*)clEnqueueMapBuffer(params.command_queue, params.positionBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(cl_double2)*params.totalParticles, 0, NULL, NULL, NULL);
	fprintf_s(positionFile, "* %f\n", currentTime);
	for (cl_uint i = 0; i < params.totalParticles; i++)
	{
		fprintf_s(positionFile, "%g %g 0 0 0 0 0 0 0 0 0\n", pointer[i].s[0], pointer[i].s[1]);
	}

	clEnqueueUnmapMemObject(params.command_queue, params.positionBuffer, pointer, NULL, NULL, NULL);
}

void outputHistogram(Parameters &params, const char* time_string)
{
	cl_double2* position = (cl_double2*)clEnqueueMapBuffer(params.command_queue, params.positionBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(cl_double2)*params.totalParticles, 0, NULL, NULL, NULL);
	
	cl_double horizLength = params.boundaries.s[1] - params.boundaries.s[0];
	cl_double vertLength = params.boundaries.s[3] - params.boundaries.s[2];
	cl_uint maxSize = (cl_uint)(sqrt(horizLength*horizLength + vertLength*vertLength)/params.binSize);
	cl_double* type1and1Histogram = new cl_double[maxSize];
	cl_double* type1and2Histogram = new cl_double[maxSize];
	cl_double* type2and2Histogram = new cl_double[maxSize];


	for (cl_uint i = 0; i < maxSize; i++)
	{
		type1and1Histogram[i] = 0;
		type1and2Histogram[i] = 0;
		type2and2Histogram[i] = 0;
	}

	cl_double2 vectors; 
	cl_double length;
	cl_uint binNumber;

	for (cl_uint i = 0; i < params.totalParticles; i++)
	{
		for (cl_uint j = i + 1; j < params.totalParticles; j++)
		{
			//This will need to be modified to accomodate moving borders
			vectors.s[0] = position[i].s[0] - position[j].s[0];// determineVector(position[i].s[0] - position[j].s[0], horizLength);
			vectors.s[1] = determineVector(position[i].s[1] - position[j].s[1], vertLength);
			length = sqrt(vectors.s[0] * vectors.s[0] + vectors.s[1] * vectors.s[1]);
			binNumber = (cl_uint)(length / params.binSize);

			if (i < params.numParticlesType1 && j < params.numParticlesType1)
			{
				type1and1Histogram[binNumber]++;
			}
			else if (i >= params.numParticlesType1 && j >= params.numParticlesType1)
			{
				type2and2Histogram[binNumber]++;
			}
			else
			{
				type1and2Histogram[binNumber]++;
			}
		}
	}

	// NEED TO DEFINE FOR EACH PARTICLE
	//double density = ((double)totalParticles) / (rightBoundary * upperBoundary);
	cl_double area;
	cl_double outerArea;
	cl_double innerArea;

	for (cl_uint i = 0; i < maxSize; i++)
	{
		outerArea = (i + 1) * params.binSize;
		outerArea *= outerArea;
		innerArea = i * params.binSize;
		innerArea *= params.binSize;

		area = 3.1415926535897932384626433832795f * (outerArea - innerArea);
		type1and1Histogram[i] /= (area);
		type1and2Histogram[i] /= (area);
		type2and2Histogram[i] /= (area);
	}

	FILE* histoFile;
	std::string fileName = time_string;
	fileName.append("/histogram_");
	fileName.append(time_string);
	fileName.append(".txt");

	fopen_s(&histoFile, fileName.c_str(), "w");

	for (cl_uint i = 0; i < maxSize; i++)
	{
		fprintf_s(histoFile, "%i %g %g %g\n", i, type1and1Histogram[i], type1and2Histogram[i], type2and2Histogram[i]);
	}

	fclose(histoFile);

	delete type1and1Histogram;
	delete type1and2Histogram;
	delete type2and2Histogram;

	clEnqueueUnmapMemObject(params.command_queue, params.positionBuffer, position, NULL, NULL, NULL);
	//With moveable boundary can't know ahead of time so have to allocate here
	/*cl_double xRange = params.boundaries.s[1] - params.boundaries.s[0];
	cl_double yRange = params.boundaries.s[3] - params.boundaries.s[2];
	cl_uint maxRange = (cl_int)(sqrt((yRange*yRange + xRange*xRange)) / (params.binSize));
	cl_int clStatus;

	cl_mem binsBuffer1to1 = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_uint)*maxRange, NULL, &clStatus);
	cl_mem binsBuffer1to2 = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_uint)*maxRange, NULL, &clStatus);
	cl_mem binsBuffer2to2 = clCreateBuffer(params.context, CL_MEM_READ_WRITE, sizeof(cl_uint)*maxRange, NULL, &clStatus);
	cl_mem adjustedBinsBufferAllTypes = clCreateBuffer(params.context, CL_MEM_WRITE_ONLY, sizeof(cl_double3)*maxRange, NULL, &clStatus);

	clSetKernelArg(params.pairCorrelationKernel, 0, sizeof(cl_mem), (void*)&binsBuffer1to1);
	clSetKernelArg(params.pairCorrelationKernel, 1, sizeof(cl_mem), (void*)&binsBuffer1to2);
	clSetKernelArg(params.pairCorrelationKernel, 2, sizeof(cl_mem), (void*)&binsBuffer2to2);
	clSetKernelArg(params.pairCorrelationKernel, 3, sizeof(cl_mem), (void*)&adjustedBinsBufferAllTypes);
	clSetKernelArg(params.pairCorrelationKernel, 4, sizeof(cl_mem), (void*)&params.positionBuffer);
	clSetKernelArg(params.pairCorrelationKernel, 5, sizeof(cl_uint), (void*)&params.totalParticles);
	clSetKernelArg(params.pairCorrelationKernel, 6, sizeof(cl_uint), (void*)&params.numParticlesType1);
	clSetKernelArg(params.pairCorrelationKernel, 7, sizeof(cl_double4), (void*)&params.boundaries);
	clSetKernelArg(params.pairCorrelationKernel, 8, sizeof(cl_uint), (void*)&maxRange);
	clSetKernelArg(params.pairCorrelationKernel, 9, sizeof(cl_double), (void*)&params.binSize);

	cout<<clEnqueueNDRangeKernel(params.command_queue, params.pairCorrelationKernel, 1, NULL, &params.globalSize, &params.localSize, 0, NULL, NULL)<<endl;

	clFlush(params.command_queue);
	clFinish(params.command_queue);

	clStatus = clReleaseMemObject(binsBuffer1to1);
	clStatus = clReleaseMemObject(binsBuffer1to2);
	clStatus = clReleaseMemObject(binsBuffer2to2);
	cl_double3* pointer = (cl_double3*)clEnqueueMapBuffer(params.command_queue, adjustedBinsBufferAllTypes, CL_TRUE, CL_MAP_READ, 0, sizeof(cl_double3)*maxRange, 0, NULL, NULL, NULL);
		
	//For some reason this won't work in the kernel, even after setting read_write
	//double area;
	//for (cl_uint i = 0; i < maxRange; i++)
	//{
		//area = 3.1415926535897932384626433832795f * (((i + 1) * params.binSize) * ((i + 1) * params.binSize) - (i * params.binSize) * (i * params.binSize));
		//pointer[i].s[0] /= (area);
		//pointer[i].s[1] /= (area);
		//pointer[i].s[2] /= (area);
	//}

	FILE* histoFile;
	string fileName = "histogram_";
	fileName.append(time_string);
	fileName.append(".txt");

	fopen_s(&histoFile, fileName.c_str(), "w");

	for (cl_uint i = 0; i < maxRange; i++)
	{
		fprintf_s(histoFile, "%i %g %g %g\n", i, pointer[i].s[0], pointer[i].s[1], pointer[i].s[2]);
	}

	fclose(histoFile);

	clEnqueueUnmapMemObject(params.command_queue, adjustedBinsBufferAllTypes, pointer, NULL, NULL, NULL);

	clStatus = clReleaseMemObject(adjustedBinsBufferAllTypes);*/
}


void outputData(Parameters &params, const char* time_string)
{	
	std::string fileName = time_string;
	fileName.append("/simStateData_");
	fileName.append(time_string);
	fileName.append(".txt");

	FILE* outputFile;
	fopen_s(&outputFile, fileName.c_str(), "w");

	cl_double2* posPointer = (cl_double2*)clEnqueueMapBuffer(params.command_queue, params.positionBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(cl_double2)*params.totalParticles, 0, NULL, NULL, NULL);
	cl_double2* oldPosPointer = (cl_double2*)clEnqueueMapBuffer(params.command_queue, params.oldPositionBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(cl_double2)*params.totalParticles, 0, NULL, NULL, NULL);

	//fprintf_s(outputFile, "%f\n", currentTime);

	for (cl_uint i = 0; i < params.totalParticles; i++)
	{
		fprintf_s(outputFile, "%f %f %f %f\n", posPointer[i].s[0], posPointer[i].s[1], oldPosPointer[i].s[0], oldPosPointer[i].s[1]);
	}
	fprintf_s(outputFile, "\n%f %f %f %f %f %f\n", params.boundaries.s[0], params.boundaries.s[1], params.oldBoundaries.s[0], params.oldBoundaries.s[1], params.origBoundaries.s[0], params.origBoundaries.s[1]);

	fclose(outputFile);

	clEnqueueUnmapMemObject(params.command_queue, params.positionBuffer, posPointer, NULL, NULL, NULL);
	clEnqueueUnmapMemObject(params.command_queue, params.oldPositionBuffer, oldPosPointer, NULL, NULL, NULL);
}

size_t get_nextpowerof2(size_t n)
{
	/*
	* Below indicates passed no is a power of 2, so return the same.
	*/
	if (!(n & (n - 1))) {
		return (n);
	}

	while (n & (n - 1)) {
		n &= n - 1;
	}
	n <<= 1;

	return n;
}

/*
void outputSettings()
{
	// Initilize Xerces.
	xercesc::XMLPlatformUtils::Initialize();

	// Pointer to our DOMImplementation.
	xercesc::DOMImplementation *pDOMImplementation = NULL;

	// Get the DOM Implementation (used for creating DOMDocuments).
	// Also see: http://www.w3.org/TR/2000/REC-DOM-Level-2-Core-20001113/core.html
	pDOMImplementation = xercesc::DOMImplementationRegistry::getDOMImplementation(xercesc::XMLString::transcode("core"));

	// Pointer to a DOMDocument.
	xercesc::DOMDocument * pDOMDocument = NULL;

	// If you're not interested in namespaces (and most of the time I'm not),
	// just use the following line instead of the one above...
	pDOMDocument = pDOMImplementation->createDocument(0, L"Hello_World", 0);
*/
	/*
	=============================================================
	Anything below this point is optional,
	although it would be a pretty boring and empty document!.
	=============================================================

	Since we are going to add nodes to the document's root,
	we need a pointer to this root node/element (which was
	previously created when the document was created).

	Note: Due to their close relationship,
	the words "node" and "element" are often used
	interchangably. Nodes are the base class,
	and Elements are the specilizations.
	*/
/*
	xercesc::DOMElement * pRootElement = NULL;
	pRootElement = pDOMDocument->getDocumentElement();


	// Create a Comment node, and then append this to the root element.
	xercesc::DOMComment * pDOMComment = NULL;
	pDOMComment = pDOMDocument->createComment(L"Dates are formatted mm/dd/yy."
		L" Don't forget XML is case-sensitive.");
	pRootElement->appendChild(pDOMComment);


	// Create an Element node, then fill in some attributes,
	// and then append this to the root element.
	xercesc::DOMElement * pDataElement = NULL;
	pDataElement = pDOMDocument->createElement(L"data");

	// Copy the current system date to a buffer, then set/create the attribute.
	wchar_t wcharBuffer[128];
	_wstrdate_s(wcharBuffer, 9);
	pDataElement->setAttribute(L"date", wcharBuffer);

	// Convert an integer to a string, then set/create the attribute.
	_itow_s(65536, wcharBuffer, 128, 10);
	pDataElement->setAttribute(L"integer", wcharBuffer);

	// Convert a doubleing-point number to a wstring,
	// then set/create the attribute.
	std::wstringstream    myStream;
	myStream.precision(8);
	myStream.setf(std::ios_base::fixed, std::ios_base::floatfield);
	myStream << 3.1415926535897932384626433832795;
	const std::wstring ws(myStream.str());
	pDataElement->setAttribute(L"double", ws.c_str());

	// Append 'pDataElement' to 'pRootElement'.
	pRootElement->appendChild(pDataElement);


	// Create an Element node, then fill in some attributes, add some text,
	// then append this to the 'pDataElement' element.
	xercesc::DOMElement * pRow = NULL;
	pRow = pDOMDocument->createElement(L"row");

	// Create some sample data.
	_itow_s(1, wcharBuffer, 128, 10);
	pRow->setAttribute(L"index", wcharBuffer);
*/
	/*
	Create a text node and append this as well. Some people
	prefer to place their data in the text node
	which is perfectly valid, others prefer to use
	the attributes. A combination of the two is quite common.
	*/
/*
	xercesc::DOMText* pTextNode = NULL;
	pTextNode = pDOMDocument->createTextNode(L"Comments and"
		L" data can also go in text nodes.");
	pRow->appendChild(pTextNode);

	pDataElement->appendChild(pRow);


	// Create another row (this time putting data and descriptions into different places).
	pRow = pDOMDocument->createElement(L"row");
	pRow->setAttribute(L"description", L"The value of PI");
	pTextNode = pDOMDocument->createTextNode(L"3.1415926535897932384626433832795");
	pRow->appendChild(pTextNode);
	pDataElement->appendChild(pRow);


	// Create another row.
	pRow = pDOMDocument->createElement(L"row");
	pRow->setAttribute(L"website", L"http://www.ted.com.au/");
	pTextNode = pDOMDocument->createTextNode(L"TED - ideas worth sharing.");
	pRow->appendChild(pTextNode);
	pDataElement->appendChild(pRow);

	pDOMDocument->release();
	xercesc::XMLPlatformUtils::Terminate();
}
*/
//! \brief Save the DOM Document to the File System.
/*
void DoOutput2File(xercesc::DOMDocument* pDOMDocument, const wchar_t * FullFilePath)
{
	xercesc::DOMImplementation *pImplement = NULL;
	xercesc::DOMLSSerializer *pSerializer = NULL; // @DOMWriter
	xercesc::LocalFileFormatTarget *pTarget = NULL;

	//Return the first registered implementation that has the desired features. In this case, we are after
	//a DOM implementation that has the LS feature... or Load/Save.

	pImplement = xercesc::DOMImplementationRegistry::getDOMImplementation(L"LS");

	//From the DOMImplementation, create a DOMWriter.
	//DOMWriters are used to serialize a DOM tree [back] into an XML document.

	pSerializer = ((xercesc::DOMImplementationLS*)pImplement)->createLSSerializer(); //@createDOMWriter();

	//This line is optional. It just sets a feature of the Serializer to make the output
	//more human-readable by inserting line-feeds, without actually inserting any new elements/nodes
	//into the DOM tree. (There are many different features to set.) Comment it out and see the difference.

	// @pSerializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, true); // 
	xercesc::DOMLSOutput *pOutput = ((xercesc::DOMImplementationLS*)pImplement)->createLSOutput();
	xercesc::DOMConfiguration *pConfiguration = pSerializer->getDomConfig();

	// Have a nice output
	if (pConfiguration->canSetParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true))
		pConfiguration->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);

	pTarget = new xercesc::LocalFileFormatTarget(FullFilePath);
	pOutput->setByteStream(pTarget);

	// @pSerializer->write(pDOMDocument->getDocumentElement(), pOutput); // missing header "<xml ...>" if used
	pSerializer->write(pDOMDocument, pOutput);

	delete pTarget;
	pOutput->release();
	pSerializer->release();
}
*/
