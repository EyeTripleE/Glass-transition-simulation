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

using namespace std;

//IF YOU'RE GOING TO USE DYNAMIC BOUNDARIES YOU MUST USE A DOUBLE INSTEAD OF AN INT

//Calculates force on each particle based on the present position of all particles
//Cutoff distance of 4.5 units.  
void calcForceMatrices(int totalParticles, double** forceBetweenX,
	double** forceBetweenY, double* positionX, double* positionY,
	int particlesType1, double& potentialEnergy, int upperBoundary,
	int rightBoundary);

//Helper method for determining force based on the type (masses) of 
//the particles interacting
char particleComboType(int particlesType1, int i, int j);

//Recalculates the acceleration of each particle based on the force
//currently applied to it.
void recalcAcceleration(double* accelerationX, double* accelerationY,
	double** forceBetweenX, double** forceBetweenY, int totalParticles,
	int numParticlesType1);

//Determines shortest X distance between particles (since we're using 
//periodic boundary condition)
double determineXVector(double positionX1, double positionX2,
	int rightBoundary);

//Determines shortest Y distance between particles (since we're using 
//periodic boundary condition)
double determineYVector(double positionY1, double positionY2,
	int upperBoundary);

//Function that performs the Euler algorithm on all particles in the set
void performEulerOperation(int totalParticles, double** forceBetweenX,
	double** forceBetweenY, double* positionX, double* positionY,
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	int upperBoundary, int rightBoundary, double* oldPositionX,
	double* oldPositionY, double* accelerationX, double* accelerationY,
	double* velocityX, double* velocityY, double currentTime,
	double timestep, ofstream& positionFile);

//Function that performs the Verlet algorithm on all particles in the set
void performVerletOperation(int totalParticles, double** forceBetweenX,
	double** forceBetweenY, double* positionX, double* positionY,
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	int upperBoundary, int rightBoundary, double* oldPositionX,
	double* oldPositionY, double* accelerationX, double* accelerationY,
	double* velocityX, double* velocityY, double currentTime,
	double timestep, ofstream& positionFile);

//Translates particles that have exited the simulation area back into the 
//Simulation region
void applyPeriodicBoundary(double &positionX, double &positionY,
	double &oldPositionX, double &oldPositionY, int rightBoundary,
	int upperBoundary);

void createHistogram(double binSize, int totalParticles, int particlesType1,
	double* positionX, double* positionY, int upperBoundary,
	int rightBoundary);

//The main function to execute the simulation
int main()
{
	clock_t tstart = clock();
	/*CONSTANTS FOR REFERENCE THESE ARE HARDCODED*/
	//constants 0.5(sigma1 + sigma2)
	//double sigma1to1 = 1; //in units of sigma
	//double sigma2to2 = 1.4; //in units of sigma
	//double sigma1to2 = 1.2; // in units of sigma
	//double epsilon = 1;// in units of epsilon
	//double massParticle1 = 1; //in units of massParticle1
	//double massParticle2 = 2; //in units of massParticle1
	//Create data files
	ofstream positionFile;
	ofstream energyFile;
	positionFile.open("position.txt");
	energyFile.open("energy.txt");

	double timestep = 0.005; //Can be arbitrarily small
	double maxTime = 5; //Can be arbitrarily long or short
	double binSize = 0.1;

	int numParticlesType1 = 1024;
	int numParticlesType2 = 0;
	int totalParticles = numParticlesType1 + numParticlesType2;

	double potentialEnergy;
	double kineticEnergy;

	//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	int rightBoundary = 50;
	int upperBoundary = 50;

	//Particle information arrays
	double *positionX = new double[totalParticles];
	double *positionY = new double[totalParticles];
	double *velocityX = new double[totalParticles];
	double *velocityY = new double[totalParticles];
	double *accelerationX = new double[totalParticles];
	double *accelerationY = new double[totalParticles];
	double *oldPositionX = new double[totalParticles];
	double *oldPositionY = new double[totalParticles];

	//Initialize to small initial velocity
	//Increase velocity by constant factor until reach target temperature.

	//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	int part = 0;
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 32; j++)
		{
			positionX[part] = i;
			positionY[part] = j;

			if (part < totalParticles / 4)
			{
				velocityX[part] = 0.01;
				velocityY[part] = 0;
			}
			else if (part < totalParticles / 2)
			{
				velocityX[part] = -0.01;
				velocityY[part] = 0;
			}
			else if (part < (totalParticles / 2 + totalParticles / 4))
			{
				velocityX[part] = 0;
				velocityY[part] = 0.01;
			}
			else
			{
				velocityX[part] = 0;
				velocityY[part] = -0.01;
			}
			part++;
		}
	}

	//Randomizing initial velocity direction
	srand((unsigned)time(NULL));
	double holderX;
	double holderY;
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
	double **forceBetweenX;
	double **forceBetweenY;
	forceBetweenX = new double*[totalParticles];
	forceBetweenY = new double*[totalParticles];
	for (int i = 0; i <= totalParticles; i++)
	{
		forceBetweenX[i] = new double[i + 1];
		forceBetweenY[i] = new double[i + 1];
	}

	//Start of simulation - setting initial time to zero
	double currentTime = 0;

	//Perform initial Euler operation to set things in motion
	positionFile << "* " << currentTime << endl;
	performEulerOperation(totalParticles, forceBetweenX, forceBetweenY,
		positionX, positionY, numParticlesType1, potentialEnergy,
		kineticEnergy, upperBoundary, rightBoundary, oldPositionX,
		oldPositionY, accelerationX, accelerationY, velocityX, velocityY,
		currentTime, timestep, positionFile);

	energyFile << currentTime << " " << kineticEnergy + potentialEnergy << " "
		<< kineticEnergy << " " << potentialEnergy << endl;

	//Main loop - performing Verlet operations for the remainder of the simulation
	for (currentTime = timestep; currentTime < maxTime; currentTime += timestep)
	{
		positionFile << "* " << currentTime << endl;
		performVerletOperation(totalParticles, forceBetweenX, forceBetweenY,
			positionX, positionY, numParticlesType1, potentialEnergy,
			kineticEnergy, upperBoundary, rightBoundary, oldPositionX,
			oldPositionY, accelerationX, accelerationY, velocityX,
			velocityY, currentTime, timestep, positionFile);

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
	cout << (double)(clock() - tstart) / CLOCKS_PER_SEC << endl;
	for (int i = 0; i < totalParticles; i++)
	{

		delete[] forceBetweenY[i];
		delete[] forceBetweenX[i];

	}
	delete[] forceBetweenY;
	delete[] forceBetweenX;

	positionFile.close();
	energyFile.close();
	return 0;
}
//End of main function

void calcForceMatrices(int totalParticles, double** forceBetweenX,
	double** forceBetweenY, double* positionX, double* positionY,
	int particlesType1, double& potentialEnergy, int upperBoundary,
	int rightBoundary)
{
	double xvector;
	double yvector;
	char pType;
	potentialEnergy = 0;
	double pythagorean;
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
				//Having the nested else's avoids calculating pType for the majority of particles.
				pType = particleComboType(particlesType1, i, j);
				if (pType == 1)
				{
					//Negative partial derivatives of the Lennard-Jones potential
					forceBetweenX[i][j] = ((24 * xvector)
						/ (pow(pythagorean, 4)))
						* ((2 / (pow(pythagorean, 3))) - 1);
					forceBetweenY[i][j] = ((24 * yvector)
						/ (pow(pythagorean, 4)))
						* ((2 / (pow(pythagorean, 3))) - 1);
					potentialEnergy += 4
						* (pow(pythagorean, -6) - pow(pythagorean, -3));
				}
				else if (pType == 2)
				{
					forceBetweenX[i][j] = ((24 * xvector * pow(1.4, 6))
						/ (pow(pythagorean, 4)))
						* ((2 * pow(1.4, 6) / (pow(pythagorean, 3))) - 1);
					forceBetweenY[i][j] = ((24 * yvector * pow(1.4, 6))
						/ (pow(pythagorean, 4)))
						* ((2 * pow(1.4, 6) / (pow(pythagorean, 3))) - 1);
					potentialEnergy += 4
						* (pow(((1.4 * 1.4) / pythagorean), 6)
						- pow(((1.4 * 1.4) / pythagorean), 3));
				}
				else
				{
					forceBetweenX[i][j] = ((24 * xvector * pow(1.2, 6))
						/ (pow(pythagorean, 4)))
						* ((2 * pow(1.2, 6) / (pow(pythagorean, 3))) - 1);
					forceBetweenY[i][j] = ((24 * yvector * pow(1.2, 6))
						/ (pow(pythagorean, 4)))
						* ((2 * pow(1.2, 6) / (pow(pythagorean, 3))) - 1);
					potentialEnergy += 4
						* (pow(((1.2 * 1.2) / pythagorean), 6)
						- pow(((1.2 * 1.2) / pythagorean), 3));
				}
			}
		}
	}
}

char particleComboType(int particlesType1, int i, int j)
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
}

void recalcAcceleration(double* accelerationX, double* accelerationY,
	double** forceBetweenX, double** forceBetweenY, int totalParticles,
	int numParticlesType1)
{
	double sumForceX;
	double sumForceY;
	for (int i = 0; i < totalParticles; i++)
	{
		sumForceX = 0;
		sumForceY = 0;
		for (int j = 0; j < totalParticles; j++)
		{
			if (j > i)
			{
				sumForceX += -1 * forceBetweenX[j][i];
				sumForceY += -1 * forceBetweenY[j][i];
			}
			else if (j < i)
			{
				sumForceX += forceBetweenX[i][j];
				sumForceY += forceBetweenY[i][j];
			}
			//Does nothing for the case j == i
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
double determineXVector(double positionX1, double positionX2, int rightBoundary)
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
		&& (positionX2 >(4.5 - (rightBoundary - positionX1))))
	{
		return (positionX1 - positionX2);
	}
	else if (positionX1 < 4.5
		&& positionX2 >(rightBoundary - (4.5 - positionX1)))
	{
		return (positionX1 - (positionX2 - rightBoundary));
	}
	else
	{
		return (positionX1 - (positionX2 + rightBoundary));
	}
}

//This function determines if a given particle should be compared to a ghost particle in the Y direction
double determineYVector(double positionY1, double positionY2, int upperBoundary)
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
		&& (positionY2 >(4.5 - (upperBoundary - positionY1))))
	{
		return (positionY1 - positionY2);
	}
	else if (positionY1 < 4.5
		&& positionY2 >(upperBoundary - (4.5 - positionY1)))
	{
		return (positionY1 - (positionY2 - upperBoundary));
	}
	else
	{
		return (positionY1 - (positionY2 + upperBoundary));
	}
}

void performEulerOperation(int totalParticles, double** forceBetweenX,
	double** forceBetweenY, double* positionX, double* positionY,
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	int upperBoundary, int rightBoundary, double* oldPositionX,
	double* oldPositionY, double* accelerationX, double* accelerationY,
	double* velocityX, double* velocityY, double currentTime,
	double timestep, ofstream& positionFile)
{
	//Calculate individual forces
	calcForceMatrices(totalParticles, forceBetweenX, forceBetweenY, positionX,
		positionY, particlesType1, potentialEnergy, upperBoundary,
		rightBoundary);

	//Sum specific forces, divide by mass to get acceleration
	recalcAcceleration(accelerationX, accelerationY, forceBetweenX,
		forceBetweenY, totalParticles, particlesType1);

	kineticEnergy = 0;

	for (int i = 0; i < totalParticles; i++)
	{
		oldPositionX[i] = positionX[i];
		oldPositionY[i] = positionY[i];

		//Euler method
		positionX[i] += (velocityX[i] * timestep);
		positionY[i] += (velocityY[i] * timestep);

		velocityX[i] += (accelerationX[i] * timestep);
		velocityY[i] += (accelerationY[i] * timestep);

		if (i < particlesType1)
			kineticEnergy += 0.5
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

void performVerletOperation(int totalParticles, double** forceBetweenX,
	double** forceBetweenY, double* positionX, double* positionY,
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	int upperBoundary, int rightBoundary, double* oldPositionX,
	double* oldPositionY, double* accelerationX, double* accelerationY,
	double* velocityX, double* velocityY, double currentTime,
	double timestep, ofstream& positionFile)
{
	//Calculate individual forces
	calcForceMatrices(totalParticles, forceBetweenX, forceBetweenY, positionX,
		positionY, particlesType1, potentialEnergy, upperBoundary,
		rightBoundary);

	//Sum specific forces, divide by mass to get acceleration
	recalcAcceleration(accelerationX, accelerationY, forceBetweenX,
		forceBetweenY, totalParticles, particlesType1);

	kineticEnergy = 0;
	double currentPosition;
	double currentDisplacement;
	double futureDisplacement;

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
			kineticEnergy += 0.5
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

void applyPeriodicBoundary(double &positionX, double &positionY,
	double &oldPositionX, double &oldPositionY, int rightBoundary,
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

void createHistogram(double binSize, int totalParticles, int particlesType1,
	double* positionX, double* positionY, int upperBoundary,
	int rightBoundary)
{
	int maxSize = (int)sqrt(
		(rightBoundary * rightBoundary + upperBoundary * upperBoundary))
		/ (binSize);
	double* type1and1Histogram = new double[maxSize];
	double* type2and2Histogram = new double[maxSize];
	double* type2and1Histogram = new double[maxSize];
	for (int i = 0; i < maxSize; i++)
	{
		type1and1Histogram[i] = 0;
		type2and1Histogram[i] = 0;
		type2and2Histogram[i] = 0;
	}

	double xVector;
	double yVector;
	double distance;
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
			binNumber = (int)(distance / binSize);

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
	double density = ((double)totalParticles)
		/ (rightBoundary * upperBoundary);
	double area;

	for (int i = 0; i < maxSize; i++)
	{
		area = 3.1415926535897932384626433832795
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
