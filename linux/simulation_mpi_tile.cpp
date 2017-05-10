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
#include <assert.h>

#include "mpi.h"

//Calculates force on each particle based on the present position of all particles
//Cutoff distance of 4.0 units. 
//Then calculates the acceleration of each particle based on the force
//currently applied to it.
//
int myxStart, myxEnd;
int myyStart, myyEnd;
int localParticles;
int col_rank, row_rank;

MPI_Comm comm_cart;
MPI_Comm comm_col;
MPI_Comm comm_row;
int dim[2], period[2], coord[2];

void calcAcceleration(double (*accelerationR)[3], double (*positionR)[3], double (*positionC)[3], 
        double totalParticles, 
	double particlesType1, double& potentialEnergy, double* boundaries);

//Function that performs the Euler algorithm on all particles in the set
void performEulerOperation(int totalParticles, double (*positionR)[3], double (*positionC)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
        double* boundaries, double (*oldPositionR)[3], double (*accelerationR)[3],
	double (*velocityR)[3], double timestep);

//Function that performs the Verlet algorithm on all particles in the set
void performVerletOperation(int totalParticles, double (*positionR)[3], double (*positionC)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
        double* boundaries, double (*oldPositionR)[3], double (*accelerationR)[3],
	double (*velocityR)[3], double timestep, double halfInvTimestep);

//Determines shortest vector from particle 1 to particle 2 (including across boundary) in one direction
double determineVector(const double p1Pos, const double p2Pos, const double size);

//Translates particles that have exited the simulation area back into the 
//Simulation region, works in one direction so do for every direction
void applyPeriodicBoundary(double &position, double &oldPosition, const double boundary);

//Outputs to position file
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[3], const int totalParticles);

void cleanup(double (*position)[3], 
	double (*oldPosition)[3], double (*velocity)[3], double (*acceleration)[3])
{
	delete[] position;
        delete[] velocity; 
        delete[] acceleration;
        delete[]  oldPosition;
}

//The main function to execute the simulation
int main(int argc, char* argv[])
{
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int reorder;

    period[0] = 0; period[1] = 0; reorder = 0;

    MPI_Dims_create(size, 2, dim);
    assert(dim[0] == dim[1]);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm_cart);
    MPI_Cart_coords(comm_cart, rank, 2, coord);

    // dim[0] is x and dim[1] is y
    MPI_Comm_split(MPI_COMM_WORLD, coord[0], rank, &comm_row);
    MPI_Comm_split(MPI_COMM_WORLD, coord[1], rank, &comm_col);

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
    if (rank == 0)
    {	
	positionFile.open("position.txt");
	energyFile.open("energy.txt");

    }

    double timestep = 0.005; //Can be arbitrarily small
    double maxTime = 1; //Can be arbitrarily long or short
    double halfInvTimestep = 0.5/timestep; //needed for Verlet

    int numParticlesType1 = atoi(argv[1]);
    int numParticlesType2 = 0;
    int totalParticles = numParticlesType1 + numParticlesType2;

    localParticles = totalParticles/dim[0]; 
    myxStart = coord[0]*localParticles;
    myxEnd   = (coord[0] + 1)*localParticles;
    myyStart = coord[1]*localParticles;
    myyEnd   = (coord[1] + 1)*localParticles;

//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
    double boundaries[] = {50.0, 50.0, 50.0}; //{upper bound of x, upper bound of y, upper bound of z}

//Particle information arrays
    double (*position)[3] = new double[totalParticles][3];
    double (*velocity)[3] = new double[totalParticles][3];
    double (*oldPosition)[3] = new double[totalParticles][3];
    double (*acceleration)[3] = new double[totalParticles][3];

    double (*positionR)[3] = new double[localParticles][3];
    double (*velocityR)[3] = new double[localParticles][3];
    double (*oldPositionR)[3] = new double[localParticles][3];
    double (*accelerationR)[3] = new double[localParticles][3];

    double (*positionC)[3] = new double[localParticles][3];
    double (*velocityC)[3] = new double[localParticles][3];
    double (*oldPositionC)[3] = new double[localParticles][3];
    double (*accelerationC)[3] = new double[localParticles][3];

    int    count = 0;

    if (rank == 0)
    {
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
				dist = determineVector(position[i][k], position[j][k], boundaries[k]);
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
		printf("Maximum iteration reached when attempting to \
                                        place particles. Try using a larger domain\n");
		cleanup(position, velocity, acceleration, oldPosition);
		exit(1);
	}
    }

    MPI_Comm_rank(comm_row, &row_rank);
    MPI_Comm_rank(comm_col, &col_rank);
 
    //------------------------------------------------------------------------- //
    // Distribute Data           						//
    MPI_Bcast(position, 3*totalParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(velocity, 3*totalParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //--------------------------------------------------------------------------// 

    for (int i = 0; i < localParticles; i++) {
        for (int j = 0; j < 3; j++) {
            velocityR[i][j] = velocity[myyStart + i][j];
            positionR[i][j] = position[myyStart + i][j];
        }
    }
    for (int i = 0; i < localParticles; i++) {
        for (int j = 0; j < 3; j++) {
            velocityC[i][j] = velocity[myxStart + i][j];
            positionC[i][j] = position[myxStart + i][j];
        }
    }

    //Start of simulation - setting initial time to zero
    double currentTime = 0;
    double totalEnergy;
    double potentialEnergy;
    double kineticEnergy;

    if (rank == 0)
    {
       outputPosition(positionFile, currentTime, position, totalParticles);
    }

	//Perform initial Euler operation to set things in motion
    performEulerOperation(totalParticles, positionR, positionC, numParticlesType1, potentialEnergy,
		kineticEnergy, boundaries, oldPositionR, accelerationR, velocityR, timestep);
    
    totalEnergy = kineticEnergy + potentialEnergy;
    if  (rank == 0)
    {
        outputPosition(positionFile, currentTime, position, totalParticles);
	energyFile << currentTime << " " << totalEnergy	<< " " << kineticEnergy << " " << potentialEnergy << "\n";
    }

    //Main loop - performing Verlet operations for the remainder of the simulation
    count = 0;
    for (currentTime = timestep; currentTime < maxTime; currentTime += timestep)
    {
        performVerletOperation(totalParticles, positionR, positionC, numParticlesType1, potentialEnergy,
	    kineticEnergy, boundaries, oldPositionR, accelerationR, velocityR, timestep, halfInvTimestep);

	count = (count + 1) % 10; //Can set print interval arbitrarily
	if(count == 0)
	{
	    totalEnergy = kineticEnergy + potentialEnergy;
            if(rank == 0)
            {
                outputPosition(positionFile, currentTime, position, totalParticles);
	        energyFile << currentTime << " " << totalEnergy	<< " " << kineticEnergy << " "\
                     << potentialEnergy << "\n";
            }
        }

        if (col_rank == row_rank)
        {
            for (int i = 0; i < localParticles; i++) {
                for (int j = 0; j < 3; j++) {
                    velocityC[i][j] = velocityR[i][j];
                    positionC[i][j] = positionR[i][j];
                }
            }
        }
        MPI_Bcast(positionC, 3*localParticles, MPI_DOUBLE, coord[1], comm_col);
    }

    if (rank == 0)
    {
	printf("%d   %g\n", size, (double)(clock() - tstart) / CLOCKS_PER_SEC);
        positionFile.close();
	energyFile.close();
    }
    
    cleanup(position, velocity, acceleration, oldPosition);

    MPI_Finalize();

    return 0;
}
//End of main function

void calcAcceleration(double (*accelerationR)[3], double (*positionR)[3], double (*positionC)[3], 
        double localParticles, double particlesType1, double& potentialEnergy, double* boundaries)
{
	double xvector, yvector, zvector;
	double sigma, sigmaPow6, sigmaPow12;
	double pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6;
	double forceX, forceY, forceZ, forceCoeff;
	potentialEnergy = 0;
	
	//Zero out acceleration, can't do this inside the mian loop because we also access the jth entry
	for(int i = 0; i < localParticles; i++)
	{
		accelerationR[i][0] = 0.0;
		accelerationR[i][1] = 0.0;
		accelerationR[i][2] = 0.0;
	}
	
        int j;

	for (int i = 0; i < localParticles; i++)
	{
	    for (j = 0; j < localParticles; j++)
	    {
            xvector = determineVector(positionR[i][0], positionC[j][0], boundaries[0]);
            yvector = determineVector(positionR[i][1], positionC[j][1], boundaries[1]);
	    zvector = determineVector(positionR[i][2], positionC[j][2], boundaries[2]);

            if (coord[0] == coord[1]  && i == j) {
		accelerationR[i][0] = 0;
		accelerationR[i][1] = 0;
		accelerationR[i][2] = 0;
                continue;
            }

	    pythagorean = ((yvector * yvector) + (xvector * xvector) + (zvector * zvector));

	    if (pythagorean < 16.0)
	    {
		//Force derived from Lennard-Jones potential
		sigma = (i < particlesType1 && j < particlesType1) ? \
                    1.0 : ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);				
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
		accelerationR[i][0] += forceX;
		accelerationR[i][1] += forceY;
		accelerationR[i][2] += forceZ;
		//Factor of four because only calculating half of interactions, so need to double PE
		potentialEnergy += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
		}
	    }
        }
}

void performEulerOperation(int totalParticles, double (*positionR)[3], double (*positionC)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
        double* boundaries, double (*oldPositionR)[3], double (*accelerationR)[3],
	double (*velocityR)[3], double timestep)
{
	calcAcceleration(accelerationR, positionR, positionC, localParticles, \
            particlesType1, potentialEnergy, boundaries);
	double dotProd;

        int llParticles = myxEnd - myxStart + 1;
        MPI_Reduce_scatter_block(MPI_IN_PLACE, accelerationR, 3*llParticles,
            MPI_DOUBLE, MPI_SUM, comm_row); 

	double tmp;
	MPI_Allreduce(&potentialEnergy, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	potentialEnergy = tmp;

	int j;
	kineticEnergy = 0;
	for (int i = 0; i < llParticles; i++)
	{
		dotProd = 0.0;
		for (j = 0; j < 3; ++j)
		{
			oldPositionR[i][j] = positionR[myxStart + i][j];
			positionR[i][j] = positionR[myxStart + i][j] + (velocityR[myxStart + i][j] * timestep);
			velocityR[i][j] = velocityR[myxStart + i][j] + (accelerationR[i][j] * timestep);
                        applyPeriodicBoundary(positionR[i][j], oldPositionR[i][j], boundaries[j]);
			dotProd += velocityR[i][j] * velocityR[i][j];
		}

		kineticEnergy += (i < particlesType1) ? 0.5 * dotProd : dotProd;
	}

        if (rank == 0)
        {
            printf("KE : %d\n", kineticEnergy);
        }

	MPI_Allreduce(&kineticEnergy, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	kineticEnergy = tmp;

	MPI_Allgather(MPI_IN_PLACE, 3*llParticles, MPI_DOUBLE, oldPositionR,
							3*llParticles, MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 3*llParticles, MPI_DOUBLE, positionR,
							3*llParticles, MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 3*llParticles, MPI_DOUBLE, velocityR,
							3*llParticles, MPI_DOUBLE, comm_row);
}

void performVerletOperation(int totalParticles, double (*positionR)[3], double (*positionC)[3],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
        double* boundaries, double (*oldPositionR)[3], double (*accelerationR)[3],
	double (*velocityR)[3], double timestep, double halfInvTimestep)
{
	calcAcceleration(accelerationR, positionR, positionC, localParticles, \
            particlesType1, potentialEnergy, boundaries);

        // Reduce scatter
        int llParticles = myxEnd - myxStart + 1;
        MPI_Reduce_scatter_block(MPI_IN_PLACE, accelerationR, 3*llParticles,
            MPI_DOUBLE, MPI_SUM, comm_row); 

	double tmp;
	MPI_Allreduce(&potentialEnergy, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	potentialEnergy = tmp;

	double currentPosition;
	double currentDisplacement;
	double futureDisplacement;
        double dtsq = timestep*timestep;
	double dotProd;
        int j;

	kineticEnergy = 0;

	for (int i = 0; i < llParticles; i++)
	{
	    // Vector Verlet Method, unroll the loops?
	    dotProd = 0;
            for(j = 0; j < 3; ++j) //Loop over all directions
	    {
		currentDisplacement = positionR[myxStart + i][j] - oldPositionR[myxStart + i][j];
		futureDisplacement = currentDisplacement + (dtsq * accelerationR[i][j]);
		currentPosition = positionR[myxStart + i][j];
		positionR[i][j] = positionR[myxStart + i][j] + futureDisplacement;
		velocityR[i][j] = (positionR[i][j] - oldPositionR[myxStart + i][j]) * (halfInvTimestep);
		oldPositionR[i][j] = currentPosition;

                applyPeriodicBoundary(positionR[i][j], oldPositionR[i][j], boundaries[j]);
				dotProd += velocityR[i][j] * velocityR[i][j];
	    }

	    kineticEnergy += (i < particlesType1) ? 0.5 * dotProd : dotProd;
	}

	MPI_Allreduce(&kineticEnergy, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	kineticEnergy = tmp;

	MPI_Allgather(MPI_IN_PLACE, 3*llParticles, MPI_DOUBLE, oldPositionR,
							3*llParticles, MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 3*llParticles, MPI_DOUBLE, positionR,
							3*llParticles, MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 3*llParticles, MPI_DOUBLE, velocityR,
							3*llParticles, MPI_DOUBLE, comm_row);
}

double determineVector(const double p1Pos, const double p2Pos, const double size)
{
    double ds = p1Pos - p2Pos;
	return (fabs(ds) > 0.5 * size) ? (ds - copysign(size, ds)) : ds;
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

//Can likely offload this to a separate thread - Otherwise would be most efficient to do in Verlet/Euler step 
//Since the values would be in the cache.
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[3], 
    const int totalParticles)
{	
	std::string str = " 0 0 0 0 0 0 0 0 0\n";
	positionFile << "* " << currentTime << "\n";
	for (int i = 0; i < totalParticles; ++i)
	{
		positionFile << position[i][0] << " " << position[i][1] << " " << position[i][2]	<< str;
	}
}
