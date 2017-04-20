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
#include <vector>
#define DIM 3

//===================================BEGIN FUNCTION HEADERS===================================
//Calculates force on each particle based on the present position of all particles
//Cutoff distance of 4.0 units. 
//Then calculates the acceleration of each particle based on the force
//currently applied to it.
void calcAcceleration(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles, 
	double particlesType1, double& potentialEnergy, double* boundaries);

//Function that performs the Euler algorithm on all particles in the set
void performEulerOperation(int totalParticles, double (*position)[DIM],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[DIM], double (*acceleration)[DIM],
	double (*velocity)[DIM], double timestep);

//Function that performs the Verlet algorithm on all particles in the set
void performVerletOperation(int totalParticles, double (*position)[DIM],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[DIM], double (*acceleration)[DIM],
	double (*velocity)[DIM], double timestep, double halfInvTimestep, double dtsq);

//Determines shortest vector from particle 1 to particle 2 (including across boundary) in one direction
double determineVectorPeriodic(const double p1Pos, const double p2Pos, const double size);

//Vector from particle 1 to particle 2 in one direction
double determineVectorFlat(const double p1Pos, const double p2Pos);

//Translates particles that have exited the simulation area back into the 
//Simulation region, works in one direction so do for every direction
void applyPeriodicBoundary(double &position, double &oldPosition, double *boundary);

//Bounce particle off wall
void applySolidBoundary(double &position, double &oldPosition, double *boundary);

//Outputs to position file
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[DIM], const int totalParticles);

//Deletes arrays and closes files
void cleanup(std::ofstream &positionFile, std::ofstream &energyFile, double (*position)[DIM], 
	double (*oldPosition)[DIM], double (*velocity)[DIM], double (*acceleration)[DIM]);

//===================================END FUNCTION HEADERS=====================================

//===================================BEGIN TREE STUFF=========================================

struct Node{
	int particleIndex;
	double mass, com[DIM];
};

class Tree{
	
	std::vector<Node> nodesArray;
	
	public:
		Tree(double (*position)[DIM], int numParticles, double boundaries[6])
		{
			std::vector<int> indices(numParticles);
			for(int i = 0; i < numParticles; i++)
				indices[i] = i;		
			buildTree(0, position, indices, boundaries);
		};

            Tree(std::vector<Node> &vec)
            {
                  nodesArray = vec;
            }

	//For parallel construction, Give each thread 1/n particles, have place into k (8 or 64) blocks. Merge the results. 
	//(Or just do that on the master thread)
	//Then spin up k processes to create
	//an array based tree. Merge the trees. 
	void buildTree(int nodeIndex, double (*position)[DIM], std::vector<int> &partIndices, double boundaries[6])
	{	
            //printf("Node index: %d, Number of particles: %ld, Boundaries: %g %g %g %g %g %g\n", nodeIndex, partIndices.size(),
            //      boundaries[0], boundaries[1], boundaries[2],boundaries[3],boundaries[4],boundaries[5]);
		if(partIndices.size() == 0)
		{
                  if(nodeIndex >= nodesArray.size())
                       nodesArray.resize(nodeIndex + 1);
			nodesArray[nodeIndex].particleIndex = -1;
		}
		else if(partIndices.size() == 1)
		{
                  if(nodeIndex >= nodesArray.size())
                       nodesArray.resize(nodeIndex + 1);
			nodesArray[nodeIndex].particleIndex = partIndices[0];
			//Assume equivalent mass, otherwise will need to modify
			nodesArray[nodeIndex].mass = 1.0;
			for(int i = 0; i < DIM; i++)
				nodesArray[nodeIndex].com[i] = position[partIndices[0]][i];		
		}
		else
		{
			//Assuming mass of one, sum mass
                  if(nodeIndex >= nodesArray.size())
                       nodesArray.resize(nodeIndex + 1);

                  nodesArray[nodeIndex].particleIndex = -2;
			nodesArray[nodeIndex].mass = partIndices.size();
			//Assuming mass of one, average location
			for(int j = 0; j < DIM; j++)
				nodesArray[nodeIndex].com[j] = 0.0;
			for(int i = 0; i < partIndices.size(); i++)
			{
				for(int j = 0; j < DIM; j++)
				{
					nodesArray[nodeIndex].com[j] += position[partIndices[i]][j]; 
				}
			}
			for(int j = 0; j < DIM; j++)
				nodesArray[nodeIndex].com[j] /= partIndices.size();
			
			std::vector<int> indicesArray[8];
			
			//Create new boundaries
			double halfX = 0.5*(boundaries[1] + boundaries[0]);
			double halfY = 0.5*(boundaries[3] + boundaries[2]);
			double halfZ = 0.5*(boundaries[5] + boundaries[4]);
			double boundaryArray[8][6] = {{boundaries[0],halfX,boundaries[2],halfY,boundaries[4],halfZ},//bottom, front, left
							      {boundaries[0],halfX,boundaries[2],halfY,halfZ,boundaries[5]},
							      {boundaries[0],halfX,halfY,boundaries[3],boundaries[4],halfZ},
							      {boundaries[0],halfX,halfY,boundaries[3],halfZ,boundaries[5]},
                                                {halfX,boundaries[1],boundaries[2],halfY,boundaries[4],halfZ},
							      {halfX,boundaries[1],boundaries[2],halfY,halfZ,boundaries[5]},
							      {halfX,boundaries[1],halfY,boundaries[3],boundaries[4],halfZ},
                                                {halfX,boundaries[1],halfY,boundaries[3],halfZ,boundaries[5]}};

			//Subdivide the indices based on boundaries
			for(int i = 0; i < partIndices.size(); i++)
			{
				if(position[partIndices[i]][0] < halfX) //Bottom
				{
					if(position[partIndices[i]][1] < halfY) //Left
					{
						if(position[partIndices[i]][2] < halfZ) //Front
						{
							indicesArray[0].push_back(partIndices[i]);
						}
						else //Back
						{
							indicesArray[1].push_back(partIndices[i]);
						}				
					}
					else //Right
					{
						if(position[partIndices[i]][2] < halfZ) //Front
						{
							indicesArray[2].push_back(partIndices[i]);
						}
						else //Back
						{
							indicesArray[3].push_back(partIndices[i]);
						}
					}
				}
				else //Top
				{
					if(position[partIndices[i]][1] < halfY) //Left
					{
						if(position[partIndices[i]][2] < halfZ) //Front
						{
							indicesArray[4].push_back(partIndices[i]);				
						}
						else //Back
						{
							indicesArray[5].push_back(partIndices[i]);
						}				
					}
					else //Right
					{
						if(position[partIndices[i]][2] < halfZ) //Front
						{
							indicesArray[6].push_back(partIndices[i]);				
						}
						else //Back
						{
							indicesArray[7].push_back(partIndices[i]);
						}
					}
				}
			}			
									
			//Build tree based on subdivisions
			for(int i = 0; i < 8; i++)
				buildTree(8*nodeIndex + (i + 1), position, indicesArray[i], boundaryArray[i]);
		}
	}
    
      void calcAcc(int nodeIndex, int partIndex, double (*position)[DIM], double acceleration[DIM], double &potentialEnergy)
      {            
            //Empty square or itself, do nothing
            if(nodesArray[nodeIndex].particleIndex != -1 && nodesArray[nodeIndex].particleIndex != partIndex)
            {
                  double vectors[DIM];
                  double pythagorean = 0;
                  for(int i = 0; i < DIM; i++)
                  {
                        vectors[i] = determineVectorFlat(position[partIndex][i], nodesArray[nodeIndex].com[i]);
                        pythagorean += vectors[i]*vectors[i];
                  }

                  double s = nodesArray[nodeIndex].mass/sqrt(pythagorean);

                  if(nodesArray[nodeIndex].particleIndex >= 0 || s < 0.5) //Calculate force
                  {
		            //Force derived from Lennard-Jones potential
                        //Guessing this should be an average of sigma values of all involved particles. Since only using type 1 just set to one here.     
		            double sigma = 1.0;//(i < particlesType1 && j < particlesType1) ? 1.0 : ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);				
		            double sigmaPow6 = sigma*sigma;//Use sigmaPow6 as to hold the square temporarily
		            sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
		            double sigmaPow12 = sigmaPow6*sigmaPow6;
		            double invPy = 1.0 / pythagorean;
		            double invPyPow3 = invPy*invPy*invPy;
		            double invPyPow4 = invPyPow3*invPy;
		            double invPyPow6 = invPyPow3*invPyPow3;
		            double forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);				

                        for(int k = 0; k < DIM; k++)
                              acceleration[k] += vectors[k]*forceCoeff;

                        potentialEnergy += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
                  }
                  else //Check children
                  {
                        for(int i = 0; i < 8; i++)
                              calcAcc(8*nodeIndex + (i + 1), partIndex, position, acceleration, potentialEnergy);
                  }   
            }   
      }
      
};

//============================END TREE STUFF============================================
//============================DECLARE TREE VERSION======================================

//Barnes-Hut version
void calcAcceleration(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles, 
	double particlesType1, double& potentialEnergy, double boundaries[6], Tree &tree);

//=======================================================================================

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
	double maxTime = 10; //Can be arbitrarily long or short
	double halfInvTimestep = 0.5/timestep; //needed for Verlet
	double dtsq = timestep*timestep;

	int numParticlesType1 = 1024;
	int numParticlesType2 = 0;
	int totalParticles = numParticlesType1 + numParticlesType2;

	double potentialEnergy, kineticEnergy, totalEnergy;

	//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	double boundaries[6] = {0.0, 50.0, 0.0, 50.0, 0.0, 50.0}; //{x_min, x_max, y_min, y_max, z_min, z_max}

	//Particle information arrays
	double (*position)[DIM] = new double[totalParticles][DIM];
	double (*velocity)[DIM] = new double[totalParticles][DIM];
	double (*oldPosition)[DIM] = new double[totalParticles][DIM];
	double (*acceleration)[DIM] = new double[totalParticles][DIM];

	std::default_random_engine generator((unsigned)time(NULL));
	std::normal_distribution<double> normalDistribution(0.0, 0.5);
      std::uniform_real_distribution<double> uniformDistribution[DIM];
      for(int i = 0; i < DIM; i++)
            uniformDistribution[i] = std::uniform_real_distribution<double>(boundaries[2*i], boundaries[2*i + 1]);
	
	for(int i = 0; i < totalParticles; i++)
            for(int j = 0; j < DIM; j++)
                  velocity[i][j] = normalDistribution(generator);

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
            for(int j = 0; j < DIM; j++)
                  position[i][j] = uniformDistribution[j](generator);
	
		for(int j = 0; j < i; j++)
		{
			norm = 0;
			for(int k = 0; k < DIM; k++)
			{	
				//dist = determineVectorPeriodic(position[i][k], position[j][k], boundaries[k]);
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
		kineticEnergy, boundaries, oldPosition, acceleration, velocity, timestep, halfInvTimestep, dtsq);

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

//O(n^2) force calculation
void calcAcceleration(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles, 
	double particlesType1, double& potentialEnergy, double boundaries[6])
{
	double xvector, yvector, zvector;
	double sigma, sigmaPow6, sigmaPow12;
	double pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6;
	double vectors[DIM], forceCoeff;
	potentialEnergy = 0;
      int i, j, k;
	
	//Zero out acceleration, can't do this inside the main loop because we also access the jth entry
	for(i = 0; i < totalParticles; i++)
            for(j = 0; j < DIM; j++)
                  acceleration[i][j] = 0.0;

	for (i = 0; i < totalParticles; i++)
	{
		for (j = 0; j < i; j++)
		{
                  pythagorean = 0;
                  for(k = 0; k < DIM; k++)
                  {
                        vectors[k] = determineVectorFlat(position[i][k], position[j][k]);
                        //vectors[k] = determineVectorPeriodic(position[i][k], position[j][k], boundaries[2*k + 1]);
                        pythagorean += vectors[k]*vectors[k];
                  }

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

                        for(k = 0; k < DIM; k++)
                        {
                              vectors[k] *= forceCoeff;
                              acceleration[i][k] += vectors[k];
                              acceleration[j][k] -= vectors[k];
                        }

				//Factor of four because only calculating half of interactions, so need to double PE
                        potentialEnergy += 4 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
			}
		}
	}

	//Because type two particles are twice as heavy
	for (i = particlesType1; i < totalParticles; ++i)
            for(j = 0; j < DIM; j++)
                  acceleration[i][j] *= 0.5;
}
/*
void calcAcceleration(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles, 
	double particlesType1, double& potentialEnergy, double boundaries[6])
{
      //Zero out current acceleration
      potentialEnergy = 0;
      //Tree tree()
	for(int i = 0; i < totalParticles; i++)
	{
            for(int j = 0; j < DIM; j++)
                  acceleration[i][j] = 0.0;		
            tree.calcAcc(0, i, position, acceleration[i], potentialEnergy);
	}
}
*/

void performEulerOperation(int totalParticles, double (*position)[DIM],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double boundaries[6], double (*oldPosition)[DIM], double (*acceleration)[DIM],
	double (*velocity)[DIM], double timestep)
{
	calcAcceleration(acceleration, position, totalParticles, particlesType1, potentialEnergy, boundaries);
	double dotProd;

	int j;
	kineticEnergy = 0;
	for (int i = 0; i < totalParticles; i++)
	{
		dotProd = 0.0;
		for (j = 0; j < DIM; ++j)
		{
			oldPosition[i][j] = position[i][j];
			position[i][j] += (velocity[i][j] * timestep);
			velocity[i][j] += (acceleration[i][j] * timestep);

                  applySolidBoundary(position[i][j], oldPosition[i][j], &boundaries[2*j]);
                  //applyPeriodicBoundary(position[i][j], oldPosition[i][j], boundaries[2*j]);
			dotProd += velocity[i][j] * velocity[i][j];
		}

		kineticEnergy += (i < particlesType1) ? 0.5 * dotProd : dotProd;
	}
}

void performVerletOperation(int totalParticles, double (*position)[DIM],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double boundaries[6], double (*oldPosition)[DIM], double (*acceleration)[DIM],
	double (*velocity)[DIM], double timestep, double halfInvTimestep, double dtsq)
{
	calcAcceleration(acceleration, position, totalParticles, particlesType1, potentialEnergy, boundaries);

	double currentPosition;
	double currentDisplacement;
	double futureDisplacement;
	double dotProd;
	int j;
	kineticEnergy = 0;

	for (int i = 0; i < totalParticles; i++)
	{
		// Vector Verlet Method, unroll the loops?
		dotProd = 0;
        	for(j = 0; j < DIM; ++j) //Loop over all directions
		{
			currentDisplacement = position[i][j] - oldPosition[i][j];
			futureDisplacement = currentDisplacement + (dtsq * acceleration[i][j]);
			currentPosition = position[i][j];
			position[i][j] = position[i][j] + futureDisplacement;
			velocity[i][j] = (position[i][j] - oldPosition[i][j]) * (halfInvTimestep);
			oldPosition[i][j] = currentPosition;

                  applySolidBoundary(position[i][j], oldPosition[i][j], &boundaries[2*j]);
                  //applyPeriodicBoundary(position[i][j], oldPosition[i][j], boundaries[2*j]);
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

//CHANGE THESE FUNCTIONS IF WANT TO SUPPORT NON-ZERO LOWER BOUND
void applyPeriodicBoundary(double &position, double &oldPosition, double *boundary)
{
	if(position < boundary[0])
	{
            double size = boundary[1] - boundary[0];
		position += size;
		oldPosition += size;
	}
  	else if(position > boundary[1])
	{
            double size = boundary[1] - boundary[0];
		position -= size;
		oldPosition -= size;
	}
}

void applySolidBoundary(double &position, double &oldPosition, double *boundary)
{
	//If pass in both axis boundaries, could do this with one if statement
	if(position < boundary[0])
	{
		position -= 2*(position - boundary[0]);
		oldPosition -= 2*(oldPosition - boundary[0]);
	}
	else if(position > boundary[1])
	{
		position -= 2*(position - boundary[1]);
		oldPosition -= 2*(oldPosition - boundary[1]);
	}
}

//Can likely offload this to a separate thread - Otherwise would be most efficient to do in Verlet/Euler step 
//Since the values would be in the cache.
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[DIM], const int totalParticles)
{	
      std::string str;
      //Fill out remaining zeros
      int i, j;
      for(i = 0; i < 3 - DIM; ++i)
            str.append("0 ");
      str.append("0 0 0 0 0 0 0 0 0\n");
	positionFile << "* " << currentTime << "\n";
	for (i = 0; i < totalParticles; ++i)
      {
            for(j = 0; j < DIM; ++j)
                  positionFile << position[i][j] << " ";
            positionFile << str;
      }
}

void cleanup(std::ofstream &positionFile, std::ofstream &energyFile, double (*position)[DIM], 
	double (*oldPosition)[DIM], double (*velocity)[DIM], double (*acceleration)[DIM])
{
	delete[] position, velocity, acceleration, oldPosition;
	positionFile.close();
	energyFile.close();
}
