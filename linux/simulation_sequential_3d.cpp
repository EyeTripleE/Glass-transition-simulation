/*
*This program simulates the movement of inert particles using Lennard-Jones potential.
*Force is calculated from the potential.  Position and velocity are updated using the Verlet
* Algorithm (with Euler algorithm initialization).
*/

#include <cmath>
#include <iostream>
#include <ctime>
#include <fstream>
#include <random>

//Calculates force on each particle based on the present position of all particles
//Cutoff distance of 4.0 units. 
//Then calculates the acceleration of each particle based on the force
//currently applied to it.
void calcAcceleration(double (*acceleration)[3], double (*position)[3], double totalParticles, 
	double particlesType1, double& potentialEnergy, double* boundaries);

//Function that performs the Euler algorithm on all particles in the set
void performEulerOperation(int totalParticles, double (*position)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[3], double (*acceleration)[3],
	double (*velocity)[3], double timestep);

//Function that performs the Verlet algorithm on all particles in the set
void performVerletOperation(int totalParticles, double (*position)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[3], double (*acceleration)[3],
	double (*velocity)[3], double timestep, double halfInvTimestep);

//Determines shortest vector from particle 1 to particle 2 (including across boundary) in one direction
double determineVectorPeriodic(const double p1Pos, const double p2Pos, const double size);

//Vector from particle 1 to particle 2 in one direction
double determineVectorFlat(const double p1Pos, const double p2Pos);

//Translates particles that have exited the simulation area back into the 
//Simulation region, works in one direction so do for every direction
void applyPeriodicBoundary(double &position, double &oldPosition, const double boundary);

//Bounce particle off wall
void applySolidBoundary(double &position, double &oldPosition, const double boundary);

//Outputs to position file
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[3], const int totalParticles);

void cleanup(std::ofstream &positionFile, std::ofstream &energyFile, double (*position)[3], 
	double (*oldPosition)[3], double (*velocity)[3], double (*acceleration)[3])
{
	delete[] position, velocity, acceleration, oldPosition;
	positionFile.close();
	energyFile.close();
}

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
	
	std::ofstream positionFile, energyFile;
	positionFile.open("position.txt");
	energyFile.open("energy.txt");

	double timestep = 0.005; //Can be arbitrarily small
	double maxTime = 1; //Can be arbitrarily long or short
	double halfInvTimestep = 0.5/timestep; //needed for Verlet

	int numParticlesType1 = 1024;
	int numParticlesType2 = 0;
	int totalParticles = numParticlesType1 + numParticlesType2;

	double potentialEnergy, kineticEnergy, totalEnergy;

	//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	double boundaries[] = {50.0, 50.0, 50.0}; //{upper bound of x, upper bound of y, upper bound of z}

	//Particle information arrays
	double (*position)[3] = new double[totalParticles][3];
	double (*velocity)[3] = new double[totalParticles][3];
	double (*oldPosition)[3] = new double[totalParticles][3];
	double (*acceleration)[3] = new double[totalParticles][3];

	std::default_random_engine generator((unsigned)time(NULL));
	std::normal_distribution<double> normalDistribution(0.0, 0.5);
	std::uniform_real_distribution<double> uniformDistributionX(0.0, boundaries[0]);
	std::uniform_real_distribution<double> uniformDistributionY(0.0, boundaries[1]);
	std::uniform_real_distribution<double> uniformDistributionZ(0.0, boundaries[2]);
	
	for(int i = 0; i < totalParticles; i++)
	{
		velocity[i][0] = normalDistribution(generator);
		velocity[i][1] = normalDistribution(generator);
		velocity[i][2] = normalDistribution(generator);		
	}

	double minSpacing = 1.0;//Be careful with this parameter, energy is really high when this is small
	double maxIter = 50;
	double norm;
	int oldi = 0;
	unsigned count = 0;
	double dist;
	//Loop for a set number of iterations or until no particles have to be moved

	for(int i = 0; i < totalParticles && count < maxIter; i++)
	{		
		count = (i != oldi) ? 0 : count + 1; //If this is a particle we've already moved, add one to count
		oldi = i;
		position[i][0] = uniformDistributionX(generator);
		position[i][1] = uniformDistributionY(generator);
		position[i][2] = uniformDistributionZ(generator);
	
		for(int j = 0; j < i; j++)
		{
			norm = 0;
			for(int k = 0; k < 3; k++)
			{	
				dist = determineVectorFlat(position[i][k], position[j][k]);
				norm += dist*dist;
			}
			norm = sqrt(norm);
			
			if(norm < minSpacing)
			{
				//Start over checking ith particle against prior particles
				i--;
				break; 
			}
		}
	}
	if(count >= maxIter)
	{
		printf("Maximum iteration reached when attempting to place particles. Try using a larger domain.\n");
		cleanup(positionFile, energyFile, position, velocity, acceleration, oldPosition);
		exit(1);
	}

	//Start of simulation - setting initial time to zero
	double currentTime = 0;

    outputPosition(positionFile, currentTime, position, totalParticles);

	//Perform initial Euler operation to set things in motion
	performEulerOperation(totalParticles, position, numParticlesType1, potentialEnergy,
	kineticEnergy, boundaries, oldPosition, acceleration, velocity, timestep);
    
	totalEnergy = kineticEnergy + potentialEnergy;
	outputPosition(positionFile, currentTime, position, totalParticles);
	energyFile << currentTime << " " << totalEnergy	<< " " << kineticEnergy << " " << potentialEnergy << "\n";

	//Main loop - performing Verlet operations for the remainder of the simulation
	count = 0;
	for (currentTime = timestep; currentTime < maxTime; currentTime += timestep)
	{
		performVerletOperation(totalParticles, position, numParticlesType1, potentialEnergy,
		kineticEnergy, boundaries, oldPosition, acceleration, velocity, timestep, halfInvTimestep);

		count = (count + 1) % 10; //Can set print interval arbitrarily
		if(count == 0)
		{
			totalEnergy = kineticEnergy + potentialEnergy;
			outputPosition(positionFile, currentTime, position, totalParticles);
			energyFile << currentTime << " " << totalEnergy << " " << kineticEnergy << " " << potentialEnergy << "\n";
		}
	}

	printf("%g\n", (double)(clock() - tstart) / CLOCKS_PER_SEC);

	cleanup(positionFile, energyFile, position, velocity, acceleration, oldPosition);
	return 0;
}
//End of main function

void calcAcceleration(double (*acceleration)[3], double (*position)[3], double totalParticles, 
	double particlesType1, double& potentialEnergy, double* boundaries)
{
	double xvector, yvector, zvector;
	double sigma, sigmaPow6, sigmaPow12;
	double pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6;
	double forceX, forceY, forceZ, forceCoeff;
	potentialEnergy = 0;
	
	//Zero out acceleration, can't do this inside the mian loop because we also access the jth entry
	for(int i = 0; i < totalParticles; i++)
	{
		acceleration[i][0] = 0.0;
		acceleration[i][1] = 0.0;
		acceleration[i][2] = 0.0;
	}
	
    int j;
	for (int i = 0; i < totalParticles; i++)
	{
		for (j = 0; j < i; j++)
		{
			xvector = determineVectorFlat(position[i][0], position[j][0]);
			yvector = determineVectorFlat(position[i][1], position[j][1]);
			zvector = determineVectorFlat(position[i][2], position[j][2]);
			pythagorean = ((yvector * yvector) + (xvector * xvector) + (zvector * zvector));

			if (pythagorean < 16.0)
			{
				//Force derived from Lennard-Jones potential
				sigma = (i < particlesType1 && j < particlesType1) ? 1.0 : ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);				
				sigmaPow6 = sigma*sigma;//Use sigmaPow6 as to hold the square temporarily
				sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
				sigmaPow12 = sigmaPow6*sigmaPow6;
				invPy = 1.0 / pythagorean;
				invPyPow3 = invPy*invPy*invPy;
				invPyPow4 = invPyPow3*invPy;
				invPyPow6 = invPyPow3*invPyPow3;
				forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);				

				forceX = forceCoeff * xvector;
				forceY = forceCoeff * yvector;
				forceZ = forceCoeff * zvector;
				acceleration[i][0] += forceX;
				acceleration[i][1] += forceY;
				acceleration[i][2] += forceZ;
				acceleration[j][0] -= forceX;
				acceleration[j][1] -= forceY;
				acceleration[j][2] -= forceZ;
				//Factor of four because only calculating half of interactions, so need to double PE
				potentialEnergy += 4 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
			}
		}
	}
	//Because type two particles are twice as heavy
	for (int i = particlesType1; i < totalParticles; ++i)
	{
		acceleration[i][0] = 0.5*acceleration[i][0];
		acceleration[i][1] = 0.5*acceleration[i][1];
		acceleration[i][2] = 0.5*acceleration[i][2];
	}   
}

void performEulerOperation(int totalParticles, double (*position)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[3], double (*acceleration)[3],
	double (*velocity)[3], double timestep)
{
	calcAcceleration(acceleration, position, totalParticles, particlesType1, potentialEnergy, boundaries);
	double dotProd;

	int j;
	kineticEnergy = 0;
	for (int i = 0; i < totalParticles; i++)
	{
		dotProd = 0.0;
		for (j = 0; j < 3; ++j)
		{
			oldPosition[i][j] = position[i][j];
			position[i][j] += (velocity[i][j] * timestep);
			velocity[i][j] += (acceleration[i][j] * timestep);

            applySolidBoundary(position[i][j], oldPosition[i][j], boundaries[j]);
			dotProd += velocity[i][j] * velocity[i][j];
		}

		kineticEnergy += (i < particlesType1) ? 0.5 * dotProd : dotProd;
	}
}

void performVerletOperation(int totalParticles, double (*position)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[3], double (*acceleration)[3],
	double (*velocity)[3], double timestep, double halfInvTimestep)
{
	calcAcceleration(acceleration, position, totalParticles, particlesType1, potentialEnergy, boundaries);

	double currentPosition;
	double currentDisplacement;
	double futureDisplacement;
	double dtsq = timestep*timestep;
	double dotProd;
	int j;
	kineticEnergy = 0;

	for (int i = 0; i < totalParticles; i++)
	{
		// Vector Verlet Method, unroll the loops?
		dotProd = 0;
        	for(j = 0; j < 3; ++j) //Loop over all directions
		{
			currentDisplacement = position[i][j] - oldPosition[i][j];
			futureDisplacement = currentDisplacement + (dtsq * acceleration[i][j]);
			currentPosition = position[i][j];
			position[i][j] = position[i][j] + futureDisplacement;
			velocity[i][j] = (position[i][j] - oldPosition[i][j]) * (halfInvTimestep);
			oldPosition[i][j] = currentPosition;

        	applySolidBoundary(position[i][j], oldPosition[i][j], boundaries[j]);
			dotProd += velocity[i][j] * velocity[i][j];
		}

		kineticEnergy += (i < particlesType1) ? 0.5 * dotProd : dotProd;
	}
}

double determineVectorPeriodic(const double p1Pos, const double p2Pos, const double size)
{
    double ds = p1Pos - p2Pos;
    return (fabs(ds) > 0.5 * size) ? (ds - copysign(size, ds)) : ds;
}

double determineVectorFlat(const double p1Pos, const double p2Pos)
{
	return p1Pos - p2Pos;
}

void applyPeriodicBoundary(double &position, double &oldPosition, const double boundary)
{
  	if (position > boundary)
	{
		position -= boundary;
		oldPosition -= boundary;
	}
	else if (position < 0.0)
	{
		position += boundary;
		oldPosition += boundary;
	}
}

void applySolidBoundary(double &position, double &oldPosition, const double boundary)
{
	//If pass in both axis boundaries, could do this with one if statement
	if(position > boundary)
	{
		position -= 2*(position - boundary);
		oldPosition -= 2*(oldPosition - boundary);
	}
	else if(position < 0.0)
	{
		position -= 2*position;
		oldPosition -= 2*oldPosition;
	}
}

//Can likely offload this to a separate thread - Otherwise would be most efficient to do in Verlet/Euler step 
//Since the values would be in the cache.
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[3], const int totalParticles)
{	
	std::string str = " 0 0 0 0 0 0 0 0 0\n";
	positionFile << "* " << currentTime << "\n";
	for (int i = 0; i < totalParticles; ++i)
	{
		positionFile << position[i][0] << " " << position[i][1] << " " << position[i][2] << str;
	}
}
