/*
 *This program simulates the movement of inert particles using Lennard-Jones potential.
 *Force is calculated from the potential.  Position and velocity are updated using the Verlet
 * Algorithm (with Euler algorithm initialization).
 */

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "KernelSourceReader.h"
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

using namespace std;

/*const char *saxpy_kernel =
"__kernel \n"
"void saxpy_kernel(float alpha, \n"
" __global float *A, \n"
" __global float *B, \n"
" __global float *C \n"
")"
"{ \n"
" //Get the index of the work-item \n"
" int index = get_global_id(0); \n"
" C[index] = alpha* A[index] + B[index]; \n"
"} \n";*/

//IF YOU'RE GOING TO USE DYNAMIC BOUNDARIES YOU MUST USE A float INSTEAD OF AN INT

//Calculates force on each particle based on the present position of all particles
//Cutoff distance of 4.5 units.  
void calcForceMatrices(int totalParticles, float* xForce,
	float* yForce, float* positionX, float* positionY,
	int particlesType1, float* potentialEnergyArray, int upperBoundary,
	int rightBoundary, cl_mem& positionXBuffer, cl_mem& positionYBuffer,
	cl_mem& forceBufferX, cl_mem& forceBufferY, cl_mem& potentialEnergyBuffer,
	cl_command_queue& command_queue, cl_kernel& forceKernel);

//Helper method for determining force based on the type (masses) of 
//the particles interacting
//char particleComboType(int particlesType1, int i, int j);

//Recalculates the acceleration of each particle based on the force
//currently applied to it.
void recalcAcceleration(float* accelerationX, float* accelerationY,
		float* xForce, float* yForce, int totalParticles,
		int numParticlesType1);

//Determines shortest X distance between particles (since we're using 
//periodic boundary condition)
float determineXVector(float positionX1, float positionX2,
	int rightBoundary);

//Determines shortest Y distance between particles (since we're using 
//periodic boundary condition)
float determineYVector(float positionY1, float positionY2,
		int upperBoundary);

//Function that performs the Euler algorithm on all particles in the set
void performEulerOperation(int totalParticles, float* xForce,
	float* yForce, float* positionX, float* positionY,
	int particlesType1, float* potentialEnergyArray, float& kineticEnergy,
	int upperBoundary, int rightBoundary, float* oldPositionX,
	float* oldPositionY, float* accelerationX, float* accelerationY,
	float* velocityX, float* velocityY, float currentTime,
	float timestep, ofstream& positionFile, cl_mem& positionBuffer, cl_mem& velocityBuffer,
	cl_mem& accelerationBuffer, cl_mem& futurePositionBuffer, cl_mem& positionXBuffer, cl_mem& positionYBuffer,
	cl_mem& forceBufferX, cl_mem& forceBufferY, cl_mem& potentialEnergyBuffer, cl_command_queue& command_queue, cl_kernel& eulerKernel, cl_kernel& forceKernel);

//Function that performs the Verlet algorithm on all particles in the set
void performVerletOperation(int totalParticles, float* xForce,
	float* yForce, float* positionX, float* positionY,
	int particlesType1, float* potentialEnergyArray, float& kineticEnergy,
	int upperBoundary, int rightBoundary, float* oldPositionX,
	float* oldPositionY, float* accelerationX, float* accelerationY,
	float* velocityX, float* velocityY, float currentTime,
	float timestep, ofstream& positionFile, cl_mem& positionBuffer, cl_mem& oldPositionBuffer,
	cl_mem& velocityBuffer, cl_mem& accelerationBuffer, cl_mem& displacementBuffer,
	cl_mem& futurePositionBuffer, cl_mem& positionXBuffer, cl_mem& positionYBuffer,
	cl_mem& forceBufferX, cl_mem& forceBufferY, cl_mem& potentialEnergyBuffer, cl_command_queue& command_queue, cl_kernel& eulerKernel, cl_kernel& forceKernel);

//Translates particles that have exited the simulation area back into the 
//Simulation region
void applyPeriodicBoundary(float &positionX, float &positionY,
		float &oldPositionX, float &oldPositionY, int rightBoundary,
		int upperBoundary);

void createHistogram(float binSize, int totalParticles, int particlesType1,
		float* positionX, float* positionY, int upperBoundary,
		int rightBoundary);

//The main function to execute the simulation
int main()
{
	/*CONSTANTS FOR REFERENCE THESE ARE HARDCODED*/
//constants 0.5(sigma1 + sigma2)
//float sigma1to1 = 1; //in units of sigma
//float sigma2to2 = 1.4; //in units of sigma
//float sigma1to2 = 1.2; // in units of sigma
//float epsilon = 1;// in units of epsilon
//float massParticle1 = 1; //in units of massParticle1
//float massParticle2 = 2; //in units of massParticle1


//Create data files
	ofstream positionFile;
	ofstream energyFile;
	positionFile.open("position.txt");
	energyFile.open("energy.txt");

	float timestep = 0.001f; //Can be arbitrarily small
	float maxTime = 0.5f; //Can be arbitrarily long or short
	float binSize = 0.1f;

	int numParticlesType1 = 512;
	int numParticlesType2 = 512;
	int totalParticles = numParticlesType1 + numParticlesType2;

	float potentialEnergy;
	float kineticEnergy;

//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	int rightBoundary = 40;
	int upperBoundary = 40;

//Particle information arrays
	float *potentialEnergyArray = new float[totalParticles];
	float *positionX = new float[totalParticles];
	float *positionY = new float[totalParticles];
	float *velocityX = new float[totalParticles];
	float *velocityY = new float[totalParticles];
	float *accelerationX = new float[totalParticles];
	float *accelerationY = new float[totalParticles];
	float *oldPositionX = new float[totalParticles];
	float *oldPositionY = new float[totalParticles];
	float *xForce = new float[totalParticles*totalParticles];
	float *yForce = new float[totalParticles*totalParticles];

	/*<<<<<<<<<<<<<<<<<<<BEGIN OPENCL STUFF>>>>>>>>>>>>>>>>>*/
	
	// Get platform and device information
	cl_platform_id* platforms = NULL;
	cl_uint num_platforms;

	//Set up the Platform
	cl_int clStatus = clGetPlatformIDs(0, NULL, &num_platforms);
	platforms = (cl_platform_id *)
	malloc(sizeof(cl_platform_id)*num_platforms);
	clStatus = clGetPlatformIDs(num_platforms, platforms, NULL);

	//Get the devices list and choose the device you want to run on
	cl_device_id *device_list = NULL;
	cl_uint num_devices;
	clStatus = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 0,
		NULL, &num_devices);
	device_list = (cl_device_id *)
	malloc(sizeof(cl_device_id)*num_devices);
	clStatus = clGetDeviceIDs(platforms[0],
	CL_DEVICE_TYPE_GPU, num_devices, device_list, NULL);

	// Create one OpenCL context for each device in the platform
	cl_context context;
	context = clCreateContext(NULL, num_devices, device_list,
		NULL, NULL, &clStatus);

	// Create a command queue
	cl_command_queue command_queue = clCreateCommandQueue(
		context, device_list[0], 0, &clStatus);
	
	//Create program from source
	cl_program programEuler = createProgram(context, device_list[0], "Euler.cl");
	cl_program programVerlet = createProgram(context, device_list[0], "Verlet.cl");
	cl_program programForce = createProgram(context, device_list[0], "ForceMatrix.cl");
	
	//Create the OpenCL kernel
	cl_kernel eulerKernel = clCreateKernel(programEuler, "euler_kernel", &clStatus);
	cl_kernel verletKernel = clCreateKernel(programVerlet, "verlet_kernel", &clStatus);
	cl_kernel forceKernel = clCreateKernel(programForce, "calc_force", &clStatus);

	//Create memory buffers on the device for each array
	cl_mem positionBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem velocityBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem accelerationBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem oldPositionBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem displacementBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem futurePositionBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem potentialEnergyBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem positionXBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem positionYBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem forceBufferX = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);
	cl_mem forceBufferY = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*totalParticles, NULL, &clStatus);

	/*<<<<<<<<<<<<<<<<<<<END OPENCL STUFF>>>>>>>>>>>>>>>>>>>>>>>>>*/
	
//Initialize to small initial velocity
//Increase velocity by constant factor until reach target temperature.

//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	int part = 0;
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 32; j++)
		{
			positionX[part] = (float)i;
			positionY[part] = (float)j;

			if (part < totalParticles / 4)
			{
				velocityX[part] = 0.01f;
				velocityY[part] = 0;
			}
			else if (part < totalParticles / 2)
			{
				velocityX[part] = -0.01f;
				velocityY[part] = 0;
			}
			else if (part < (totalParticles / 2 + totalParticles / 4))
			{
				velocityX[part] = 0;
				velocityY[part] = 0.01f;
			}
			else
			{
				velocityX[part] = 0;
				velocityY[part] = -0.01f;
			}
			part++;
		}
	}

	//Randomizing initial velocity direction
	srand((unsigned) time(NULL));
	float holderX;
	float holderY;
	int locationY;
	int locationX;
	for (int i = 0; i < totalParticles; i++)
	{
		//Select two random particles in array
		locationX = rand() % totalParticles;
		locationY = rand() % totalParticles;
		//Swap velocity
		holderX = velocityX[i];
		holderY = velocityY[i];
		velocityX[i] = velocityX[locationX];
		velocityY[i] = velocityY[locationY];
		velocityX[locationX] = holderX;
		velocityY[locationY] = holderY;

		//Select single random particle in array
		locationX = rand() % totalParticles;
		//Swap position
		holderX = positionX[i];
		holderY = positionY[i];
		positionX[i] = positionX[locationX];
		positionY[i] = positionY[locationX];
		positionX[locationX] = holderX;
		positionY[locationX] = holderY;
	}

//Setting up force matrices
	/*float **forceBetweenX;
	float **forceBetweenY;
	forceBetweenX = new float*[totalParticles];
	forceBetweenY = new float*[totalParticles];
	for (int i = 0; i < totalParticles; i++)
	{
		forceBetweenX[i] = new float[i + 1];
		forceBetweenY[i] = new float[i + 1];
	}*/

//Start of simulation - setting initial time to zero
	float currentTime = 0;

//Perform initial Euler operation to set things in motion
	positionFile << "* " << currentTime << endl;

	performEulerOperation(totalParticles, xForce,
		yForce, positionX, positionY,
		numParticlesType1, potentialEnergyArray, kineticEnergy,
		upperBoundary, rightBoundary, oldPositionX,
		oldPositionY, accelerationX, accelerationY,
		velocityX, velocityY, currentTime,
		timestep, positionFile, positionBuffer, velocityBuffer,
		accelerationBuffer, futurePositionBuffer, positionXBuffer, positionYBuffer,
		forceBufferX, forceBufferY, potentialEnergyBuffer, command_queue,
		eulerKernel, forceKernel);

	//performsaxby(timestep, positionX, velocityX, accelerationX, command_queue, kernel, positionBuffer, velocityBuffer, accelerationBuffer);

	potentialEnergy = 0;
	for (int i = 0; i < totalParticles; i++)
	{
		potentialEnergy += potentialEnergyArray[i];
	}

	energyFile << currentTime << " " << kineticEnergy + potentialEnergy << " "
			<< kineticEnergy << " " << potentialEnergy << endl;

	//Main loop - performing Verlet operations for the remainder of the simulation
	for (currentTime = timestep; currentTime < maxTime; currentTime += timestep)
	{
		positionFile << "* " << currentTime << endl;
		performVerletOperation(totalParticles, xForce,
			yForce, positionX, positionY,
			numParticlesType1, potentialEnergyArray, kineticEnergy,
			upperBoundary, rightBoundary, oldPositionX,
			oldPositionY, accelerationX, accelerationY,
			velocityX, velocityY, currentTime,
			timestep, positionFile, positionBuffer, oldPositionBuffer,
			velocityBuffer, accelerationBuffer, displacementBuffer,
			futurePositionBuffer, positionXBuffer, positionYBuffer,
			forceBufferX, forceBufferY, potentialEnergyBuffer,
			command_queue, verletKernel, forceKernel);

		potentialEnergy = 0;
		for (int i = 0; i < totalParticles; i++)
		{
			potentialEnergy += potentialEnergyArray[i];
		}
	
		energyFile << currentTime << " " << kineticEnergy + potentialEnergy
				<< " " << kineticEnergy << " " << potentialEnergy << endl;
	}

	createHistogram(binSize, totalParticles, numParticlesType1, positionX,
			positionY, upperBoundary, rightBoundary);

	//Deleting arrays					
	delete positionX;
	delete positionY;
	delete velocityX;
	delete velocityY;
	delete accelerationX;
	delete accelerationY;
	delete oldPositionX;
	delete oldPositionY;
	delete xForce;
	delete yForce;
	delete potentialEnergyArray;
	/*for (int i = 0; i < totalParticles; i++)
	{
		delete[] forceBetweenY[i];
		delete[] forceBetweenX[i];
	}
	delete[] forceBetweenY;
	delete[] forceBetweenX;*/

	clStatus = clReleaseKernel(eulerKernel);
	clStatus = clReleaseProgram(programEuler);
	clStatus = clReleaseKernel(verletKernel);
	clStatus = clReleaseProgram(programVerlet);
	clStatus = clReleaseKernel(forceKernel);
	clStatus = clReleaseProgram(programForce);
	clStatus = clReleaseMemObject(positionBuffer);
	clStatus = clReleaseMemObject(oldPositionBuffer);
	clStatus = clReleaseMemObject(velocityBuffer);
	clStatus = clReleaseMemObject(accelerationBuffer);
	clStatus = clReleaseMemObject(displacementBuffer);
	clStatus = clReleaseMemObject(futurePositionBuffer);
	clStatus = clReleaseCommandQueue(command_queue);
	clStatus = clReleaseContext(context);
	free(platforms);
	free(device_list);

	positionFile.close();
	energyFile.close();
}
//End of main function

void calcForceMatrices(int totalParticles, float* xForce,
		float* yForce, float* positionX, float* positionY,
		int particlesType1, float* potentialEnergyArray, int upperBoundary,
		int rightBoundary, cl_mem& positionXBuffer, cl_mem& positionYBuffer,
		cl_mem& forceBufferX, cl_mem& forceBufferY, cl_mem& potentialEnergyBuffer,
		cl_command_queue& command_queue, cl_kernel& forceKernel) 
{
	
	//Enqueue data
	clEnqueueWriteBuffer(command_queue, positionXBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionX, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, positionYBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionY, 0, NULL, NULL);

	//Set arguments for kernel
	clSetKernelArg(forceKernel, 0, sizeof(cl_mem), (void *)&positionXBuffer);
	clSetKernelArg(forceKernel, 1, sizeof(cl_mem), (void *)&positionYBuffer);
	clSetKernelArg(forceKernel, 2, sizeof(cl_mem), (void *)&xForce);
	clSetKernelArg(forceKernel, 3, sizeof(cl_mem), (void *)&yForce);
	clSetKernelArg(forceKernel, 4, sizeof(cl_mem), (void *)&potentialEnergyBuffer);
	clSetKernelArg(forceKernel, 5, sizeof(int), (void *)&rightBoundary);
	clSetKernelArg(forceKernel, 6, sizeof(int), (void *)&upperBoundary);
	clSetKernelArg(forceKernel, 7, sizeof(int), (void *)&totalParticles);
	clSetKernelArg(forceKernel, 8, sizeof(int), (void *)&particlesType1);

	size_t global_size = totalParticles;
	size_t local_size = 64;
	clEnqueueNDRangeKernel(command_queue, forceKernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);

	//Read data back
	clEnqueueReadBuffer(command_queue, potentialEnergyBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, potentialEnergyArray, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, forceBufferX, CL_TRUE, 0, sizeof(float)*totalParticles, xForce, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, forceBufferY, CL_TRUE, 0, sizeof(float)*totalParticles, yForce, 0, NULL, NULL);

	clFlush(command_queue);
	clFinish(command_queue);



	/*float xvector;
	float yvector;
	potentialEnergy = 0;
	float pythagorean;
	for (int i = 0; i < totalParticles; i++)
	{
		for (int j = 0; j < i; j++)
		{
			xvector = determineXVector(positionX[i], positionX[j],
					rightBoundary);
			yvector = determineYVector(positionY[i], positionY[j],
					upperBoundary);
			pythagorean = ((yvector * yvector) + (xvector * xvector));

			if (sqrt(pythagorean) > 4.5)
			{
				forceBetweenX[i][j] = 0;
				forceBetweenY[i][j] = 0;
			}
			else
			{
				char pType = particleComboType(particlesType1, i, j);
			
				float sigma;
				if (pType == 1)
				{
					sigma = 1;
				}
				else if (pType == 2)
				{
					sigma = 1.4f;
				}
				else
				{
					sigma = 1.2f;
				}

				forceBetweenX[i][j] = ((24 * xvector * pow(sigma, 6))
					/ (pow(pythagorean, 4)))
					* ((2 * pow(sigma, 6) / (pow(pythagorean, 3))) - 1);
				forceBetweenY[i][j] = ((24 * yvector * pow(sigma, 6))
					/ (pow(pythagorean, 4)))
					* ((2 * pow(sigma, 6) / (pow(pythagorean, 3))) - 1);
				potentialEnergy += 4
					* (pow(((sigma * sigma) / pythagorean), 6)
					- pow(((sigma * sigma) / pythagorean), 3));
			}
		}
	}*/
}

/*char particleComboType(int particlesType1, int i, int j)
{
	if (i < particlesType1 && j < particlesType1)
	{
		return 1;
	}
	else if (i >= particlesType1 && j >= particlesType1)
	{
		return 2;
	}
	else
	{
		return 12;
	}
}*/

void recalcAcceleration(float* accelerationX, float* accelerationY,
		float* xForce, float* yForce, int totalParticles,
		int numParticlesType1)
{
	float sumForceX;
	float sumForceY;
	int location;
	for (int i = 0; i < totalParticles; i++)
	{
		sumForceX = 0;
		sumForceY = 0;
		for (int j = 0; j < totalParticles; j++)
		{
			if (i != j)
			{
				location = i*totalParticles + j;
				sumForceX += xForce[location];
				sumForceY += yForce[location];
			}
			/*if (j > i)
			{
				sumForceX += -1 * xForce[j][i];
				sumForceY += -1 * yForce[j][i];
			}
			else if (j < i)
			{
				sumForceX += xForce[i][j];
				sumForceY += yForce[i][j];
			}
			//Does nothing for the case j == i*/
		}
		accelerationX[i] = sumForceX;
		accelerationY[i] = sumForceY;
	}

	//Because type two particles are twice as heavy
	for (int i = numParticlesType1; i < totalParticles; i++)
	{
		accelerationX[i] = accelerationX[i] / 2;
		accelerationY[i] = accelerationY[i] / 2;
	}
}


//Determines if a given particle should be compared to a ghost particle in the X direction
float determineXVector(float positionX1, float positionX2, int rightBoundary)
{
	//This can be simplified
	//Changed || to &&
	if ((positionX1 >= 4.5) && (positionX1 <= rightBoundary - 4.5))
	{
		return (positionX1 - positionX2);
	}
	else if (positionX1 < 4.5
			&& positionX2 < (rightBoundary - (4.5 - positionX1)))
	{
		return (positionX1 - positionX2);
	}
	else if ((positionX1 > rightBoundary - 4.5)
			&& (positionX2 > (4.5 - (rightBoundary - positionX1))))
	{
		return (positionX1 - positionX2);
	}
	else if (positionX1 < 4.5
			&& positionX2 > (rightBoundary - (4.5 - positionX1)))
	{
		return (positionX1 - (positionX2 - rightBoundary));
	}
	else
	{
		return (positionX1 - (positionX2 + rightBoundary));
	}
}

//This function determines if a given particle should be compared to a ghost particle in the Y direction
float determineYVector(float positionY1, float positionY2, int upperBoundary)
{
	//This can be simplified
	//Changed || to &&
	if ((positionY1 >= 4.5) && (positionY1 <= upperBoundary - 4.5))
	{
		return (positionY1 - positionY2);
	}
	else if (positionY1 < 4.5
			&& positionY2 < (upperBoundary - (4.5 - positionY1)))
	{
		return (positionY1 - positionY2);
	}
	else if ((positionY1 > upperBoundary - 4.5)
			&& (positionY2 > (4.5 - (upperBoundary - positionY1))))
	{
		return (positionY1 - positionY2);
	}
	else if (positionY1 < 4.5
			&& positionY2 > (upperBoundary - (4.5 - positionY1)))
	{
		return (positionY1 - (positionY2 - upperBoundary));
	}
	else
	{
		return (positionY1 - (positionY2 + upperBoundary));
	}
}

void performEulerOperation(int totalParticles, float* xForce,
		float* yForce, float* positionX, float* positionY,
		int particlesType1, float* potentialEnergyArray, float& kineticEnergy,
		int upperBoundary, int rightBoundary, float* oldPositionX,
		float* oldPositionY, float* accelerationX, float* accelerationY,
		float* velocityX, float* velocityY, float currentTime,
		float timestep, ofstream& positionFile, cl_mem& positionBuffer, cl_mem& velocityBuffer,
		cl_mem& accelerationBuffer, cl_mem& futurePositionBuffer, cl_mem& positionXBuffer, cl_mem& positionYBuffer,
		cl_mem& forceBufferX, cl_mem& forceBufferY, cl_mem& potentialEnergyBuffer, cl_command_queue& command_queue,
		cl_kernel& eulerKernel, cl_kernel& forceKernel)
{
//Calculate individual forces
	calcForceMatrices(totalParticles, xForce,
		yForce, positionX, positionY,
		particlesType1, potentialEnergyArray, upperBoundary,
		rightBoundary, positionXBuffer, positionYBuffer,
		forceBufferX, forceBufferY, potentialEnergyBuffer,
		command_queue, forceKernel);

//Sum specific forces, divide by mass to get acceleration
	recalcAcceleration(accelerationX, accelerationY, xForce,
			yForce, totalParticles, particlesType1);

	kineticEnergy = 0;

	//Start with the X direction

	//Enqueue data
	clEnqueueWriteBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionX, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, velocityBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, velocityX, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, accelerationBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, accelerationX, 0, NULL, NULL);

	//Set arguments for kernel
	clSetKernelArg(eulerKernel, 0, sizeof(float), (void *)&timestep);
	clSetKernelArg(eulerKernel, 1, sizeof(cl_mem), (void *)&positionBuffer);
	clSetKernelArg(eulerKernel, 2, sizeof(cl_mem), (void *)&velocityBuffer);
	clSetKernelArg(eulerKernel, 3, sizeof(cl_mem), (void *)&accelerationBuffer);
	clSetKernelArg(eulerKernel, 4, sizeof(cl_mem), (void *)&futurePositionBuffer);

	size_t global_size = totalParticles;
	size_t local_size = 64;
	clEnqueueNDRangeKernel(command_queue, eulerKernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);

	//Read data back
	clEnqueueReadBuffer(command_queue, futurePositionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionX, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, oldPositionX, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, velocityBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, velocityX, 0, NULL, NULL);

	clFlush(command_queue);
	clFinish(command_queue);

	//Now for the Y direction

	//Enqueue data
	clEnqueueWriteBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionY, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, velocityBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, velocityY, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, accelerationBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, accelerationY, 0, NULL, NULL);

	//Set arguments for kernel
	clSetKernelArg(eulerKernel, 0, sizeof(float), (void *)&timestep);
	clSetKernelArg(eulerKernel, 1, sizeof(cl_mem), (void *)&positionBuffer);
	clSetKernelArg(eulerKernel, 2, sizeof(cl_mem), (void *)&velocityBuffer);
	clSetKernelArg(eulerKernel, 3, sizeof(cl_mem), (void *)&accelerationBuffer);
	clSetKernelArg(eulerKernel, 4, sizeof(cl_mem), (void *)&futurePositionBuffer);

	clEnqueueNDRangeKernel(command_queue, eulerKernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);

	//Read data back
	clEnqueueReadBuffer(command_queue, futurePositionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionY, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, oldPositionY, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, velocityBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, velocityY, 0, NULL, NULL);

	clFlush(command_queue);
	clFinish(command_queue);

	/*for (int i = 0; i < totalParticles; i++)
	{
		cout << i << " " << positionX[i] << " " << oldPositionX[i] << " " << velocityX[i] << endl;
		cout << i << " " << positionY[i] << " " << oldPositionY[i] << " " << velocityY[i] << endl;
	}*/

	for (int i = 0; i < totalParticles; i++)
	{
		/*oldPositionX[i] = positionX[i];
		oldPositionY[i] = positionY[i];

		//Euler method
		positionX[i] += (velocityX[i] * timestep);
		positionY[i] += (velocityY[i] * timestep);

		velocityX[i] += (accelerationX[i] * timestep);
		velocityY[i] += (accelerationY[i] * timestep);*/

		if (i < particlesType1)
			kineticEnergy += 0.5f
					* ((velocityX[i] * velocityX[i])
							+ (velocityY[i] * velocityY[i]));
		else
			kineticEnergy += ((velocityX[i] * velocityX[i])
					+ (velocityY[i] * velocityY[i]));

		applyPeriodicBoundary(positionX[i], positionY[i], oldPositionX[i],
				oldPositionY[i], rightBoundary, upperBoundary);
		positionFile << positionX[i] << " " << positionY[i]
				<< " 0 0 0 0 0 0 0 0 0 0" << '\n';
	}
}

void performVerletOperation(int totalParticles, float* xForce,
		float* yForce, float* positionX, float* positionY,
		int particlesType1, float* potentialEnergyArray, float& kineticEnergy,
		int upperBoundary, int rightBoundary, float* oldPositionX,
		float* oldPositionY, float* accelerationX, float* accelerationY,
		float* velocityX, float* velocityY, float currentTime,
		float timestep, ofstream& positionFile, cl_mem& positionBuffer, cl_mem& oldPositionBuffer, 
		cl_mem& velocityBuffer, cl_mem& accelerationBuffer, cl_mem& displacementBuffer,
		cl_mem& futurePositionBuffer, cl_mem& positionXBuffer, cl_mem& positionYBuffer,
		cl_mem& forceBufferX, cl_mem& forceBufferY, cl_mem& potentialEnergyBuffer, 
		cl_command_queue& command_queue, cl_kernel& verletKernel, cl_kernel& forceKernel)
{
	//Calculate individual forces
	calcForceMatrices(totalParticles, xForce,
		yForce, positionX, positionY,
		particlesType1, potentialEnergyArray, upperBoundary,
		rightBoundary, positionXBuffer, positionYBuffer,
		forceBufferX, forceBufferY, potentialEnergyBuffer,
		command_queue, forceKernel);

//Sum specific forces, divide by mass to get acceleration
	recalcAcceleration(accelerationX, accelerationY, xForce,
			yForce, totalParticles, particlesType1);

	kineticEnergy = 0;
	
	//Start with X direction

	//Write information to CL buffers
	/*clEnqueueWriteBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionX, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, oldPositionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, oldPositionX, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, accelerationBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, accelerationX, 0, NULL, NULL);
	
	//Set arguments for kernel
	clSetKernelArg(kernel, 0, sizeof(float), (void *) &timestep);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &positionBuffer);
	clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &oldPositionBuffer);	
	clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &velocityBuffer);
	clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &accelerationBuffer);
	clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &displacementBuffer);
	clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &futurePositionBuffer);


	//Execute the Kernel
	size_t global_size = totalParticles;
	size_t local_size = 64;
	clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);

	//read data back into arrays
	//future position becomes current position, current position becomes old position

	cout << positionX[1] << endl;

	clEnqueueReadBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, oldPositionX, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, futurePositionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionX, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, velocityBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, velocityX, 0, NULL, NULL);

	cout << positionX[1] << endl;

	clFlush(command_queue);
	clFinish(command_queue);

	//Now do for the Y direction*/

	//Write information to CL buffers
	/*clEnqueueWriteBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionY, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, oldPositionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, oldPositionY, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, accelerationBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, accelerationY, 0, NULL, NULL);

	//Set arguments for kernel
	clSetKernelArg(kernel, 0, sizeof(float), (void *)& timestep);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)& positionBuffer);
	clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)& oldPositionBuffer);
	clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)& velocityBuffer);
	clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)& accelerationBuffer);
	clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)& displacementBuffer);
	clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)& futurePositionBuffer);

	//Execute the Kernel
	clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);

	//read data back into arrays
	clEnqueueReadBuffer(command_queue, positionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, oldPositionY, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, futurePositionBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, positionY, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, velocityBuffer, CL_TRUE, 0, sizeof(float)*totalParticles, velocityY, 0, NULL, NULL);

	clFlush(command_queue);
	clFinish(command_queue);*/

	/*for (int i = 0; i < totalParticles; i++)
	{
		cout << i << " " << positionX[i] << " " << oldPositionX[i] << " " << velocityX[i] << endl;
		cout << i << " " << positionY[i] << " " << oldPositionY[i] << " " << velocityY[i] << endl;
	}*/


	float currentDisplacement;
	float currentPosition;
	float futureDisplacement;
	for (int i = 0; i < totalParticles; i++)
	{
		//Leapfrog Verlet method - unstable
		/*velocityX[i] = velocityX[i] + accelerationX[i]*timestep;
		 positionX[i] = positionX[i] + timestep*velocityX[i];
		 velocityY[i] = velocityY[i] + accelerationY[i]*timestep;
		 positionY[i] = positionY[i] + timestep*velocityY[i];*/

		// Vector Verlet Method
		currentDisplacement = positionX[i] - oldPositionX[i];
		futureDisplacement = currentDisplacement
				+ (timestep * timestep * accelerationX[i]);
		currentPosition = positionX[i];
		positionX[i] = positionX[i] + futureDisplacement;
		velocityX[i] = (positionX[i] - oldPositionX[i]) / (2 * timestep);
		oldPositionX[i] = currentPosition;

		currentDisplacement = positionY[i] - oldPositionY[i];
		futureDisplacement = currentDisplacement
				+ (timestep * timestep * accelerationY[i]);
		currentPosition = positionY[i];
		positionY[i] = positionY[i] + futureDisplacement;
		velocityY[i] = (positionY[i] - oldPositionY[i]) / (2 * timestep);
		oldPositionY[i] = currentPosition;

		//Verlet method
		/*currentPosition = positionX[i];
		 positionX[i] = 2*(currentPosition) - oldPositionX[i] + (timestep*timestep*accelerationX[i]);
		 velocityX[i] = (positionX[i] - oldPositionX[i])/(2*timestep);
		 oldPositionX[i] = currentPosition;

		 currentPosition = positionY[i];
		 positionY[i] = 2*(currentPosition) - oldPositionY[i] + (timestep*timestep*accelerationY[i]);
		 velocityY[i] = (positionY[i] - oldPositionY[i])/(2*timestep);
		 oldPositionY[i] = currentPosition;*/

		if (i < particlesType1)
			kineticEnergy += 0.5f
					* ((velocityX[i] * velocityX[i])
							+ (velocityY[i] * velocityY[i]));
		else
			kineticEnergy += ((velocityX[i] * velocityX[i])
					+ (velocityY[i] * velocityY[i]));

		//cout << kineticEnergy << endl;

		applyPeriodicBoundary(positionX[i], positionY[i], oldPositionX[i],
				oldPositionY[i], rightBoundary, upperBoundary);
		positionFile << positionX[i] << " " << positionY[i]
				<< " 0 0 0 0 0 0 0 0 0 0" << '\n';
	}
}

void applyPeriodicBoundary(float &positionX, float &positionY,
		float &oldPositionX, float &oldPositionY, int rightBoundary,
		int upperBoundary)
{
	while (positionX > rightBoundary)
	{
		positionX -= rightBoundary;
		oldPositionX -= rightBoundary;
	}
	while (positionX < 0)
	{
		positionX += rightBoundary;
		oldPositionX += rightBoundary;
	}
	while (positionY > upperBoundary)
	{
		positionY -= upperBoundary;
		oldPositionY -= upperBoundary;
	}
	while (positionY < 0)
	{
		positionY += upperBoundary;
		oldPositionY += upperBoundary;
	}
}

void createHistogram(float binSize, int totalParticles, int particlesType1,
		float* positionX, float* positionY, int upperBoundary,
		int rightBoundary)
{
	int maxSize = (int) (sqrt(
			((float)(rightBoundary * rightBoundary + upperBoundary * upperBoundary)))
			/ (binSize));
	float* type1and1Histogram = new float[maxSize];
	float* type2and2Histogram = new float[maxSize];
	float* type2and1Histogram = new float[maxSize];
	for (int i = 0; i < maxSize; i++)
	{
		type1and1Histogram[i] = 0;
		type2and1Histogram[i] = 0;
		type2and2Histogram[i] = 0;
	}

	float xVector;
	float yVector;
	float distance;
	int binNumber;

	for (int i = 0; i < totalParticles; i++)
	{
		for (int j = i + 1; j < totalParticles; j++)
		{
			xVector = determineXVector(positionX[i], positionX[j],
					rightBoundary);
			yVector = determineYVector(positionY[i], positionY[j],
					upperBoundary);
			distance = sqrt(xVector * xVector + yVector * yVector);
			binNumber = (int) (distance / binSize);

			if (i < particlesType1 && j < particlesType1)
			{
				type1and1Histogram[binNumber] += 1;
			}
			else if (i >= particlesType1 && j >= particlesType1)
			{
				type2and2Histogram[binNumber] += 1;
			}
			else
			{
				type2and1Histogram[binNumber] += 1;
			}
		}
	}

	// NEED TO DEFINE FOR EACH PARTICLE
	float density = ((float) totalParticles)
			/ (rightBoundary * upperBoundary);
	float area;

	for (int i = 0; i < maxSize; i++)
	{
		area = 3.1415926535897932384626433832795f
				* (((i + 1) * binSize) * ((i + 1) * binSize)
						- (i * binSize) * (i * binSize));
		type1and1Histogram[i] *= (1 / area);
		type2and1Histogram[i] *= (1 / area);
		type2and2Histogram[i] *= (1 / area);
	}

	ofstream histoFile;
	histoFile.open("histoFile.txt");
	for (int i = 0; i < maxSize; i++)
	{
		histoFile << i << " " << type1and1Histogram[i] << " "
				<< type2and1Histogram[i] << " " << type2and2Histogram[i]
				<< endl;
	}
	histoFile.close();

	delete type1and1Histogram;
	delete type2and1Histogram;
	delete type2and2Histogram;
}