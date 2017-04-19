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

//===================================BEGIN TREE STUFF=========================================

struct Node{
	int particleIndex;
	double mass, com[DIM];
};

class Tree{
	
	std::vector<Node> nodesArray;
	
	public:
		Tree(double (*position)[DIM], int numParticles, double boundaries[DIM])
		{
			std::vector<int> indices(numParticles);
			for(int i = 0; i < numParticles; i++)
				indices[i] = i;
			double bound[6] = {0.0, boundaries[0], 0.0, boundaries[1], 0.0, boundaries[2]};		
			buildTree(0, position, indices, bound);
		};
	
	private:
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


      /*      
      void calcAcc(int nodeIndex, int partIndex, double (*position)[DIM], double acceleration[DIM])
      {
            
            
            double distance = 0;
            double val;
            for(int i = 0; i < DIM; i++)
                  val = determineVectorFlat(position[partIndex][i], nodesArray[nodeIndex].com[i]) 
                  distance += val*val;
            distance = sqrt(distance);

            double s = nodesArray[nodeIndex].mass/distance;

            if(nodesArray[nodeIndex].particleIndex == -1)
            {
                  return;
            }
            else if(nodesArray[nodeIndex].particleIndex >= 0)
            {
                  double xvector = determineVectorFlat(position[i][0], position[j][0]);
			double yvector = determineVectorFlat(position[i][1], position[j][1]);
			double zvector = determineVectorFlat(position[i][2], position[j][2]);
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
				//Factor of four because only calculating half of interactions, so need to double PE
				potentialEnergy += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
			}
            }
            else
            {
                  if()
            }

            
            
            //mass/distance < 0.5 then ignore internal structure
            //Check if ratio is less than cutoff
            //If it is not -> look at children unless we've reached a leaf node (in which case calculate), if at leaf node make sure not looking at self
            //If it is greater, calculate the accleration
      }
      */
      
};

//============================END TREE STUFF============================================

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
	double (*velocity)[DIM], double timestep, double halfInvTimestep);

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
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[DIM], const int totalParticles);

void cleanup(std::ofstream &positionFile, std::ofstream &energyFile, double (*position)[DIM], 
	double (*oldPosition)[DIM], double (*velocity)[DIM], double (*acceleration)[DIM])
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
	double maxTime = 10; //Can be arbitrarily long or short
	double halfInvTimestep = 0.5/timestep; //needed for Verlet

	int numParticlesType1 = 1024;
	int numParticlesType2 = 0;
	int totalParticles = numParticlesType1 + numParticlesType2;

	double potentialEnergy, kineticEnergy, totalEnergy;

	//THIS MUST BE MODIFIED IF YOU USE MORE PARTICLES
	double boundaries[] = {50.0, 50.0, 50.0}; //{upper bound of x, upper bound of y, upper bound of z}

	//Particle information arrays
	double (*position)[DIM] = new double[totalParticles][DIM];
	double (*velocity)[DIM] = new double[totalParticles][DIM];
	double (*oldPosition)[DIM] = new double[totalParticles][DIM];
	double (*acceleration)[DIM] = new double[totalParticles][DIM];

	std::default_random_engine generator((unsigned)time(NULL));
	std::normal_distribution<double> normalDistribution(0.0, 0.5);
      std::uniform_real_distribution<double> uniformDistribution[DIM];
      for(int i = 0; i < DIM; i++)
            uniformDistribution[i] = std::uniform_real_distribution<double>(0.0, boundaries[i]);
	
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

      //Tree tree(position, totalParticles, boundaries);
      
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

void calcAcceleration(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles, 
	double particlesType1, double& potentialEnergy, double* boundaries)
{
	double xvector, yvector, zvector;
	double sigma, sigmaPow6, sigmaPow12;
	double pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6;
	double vectors[DIM], forceCoeff;
	potentialEnergy = 0;
      int j, k;
	
	//Zero out acceleration, can't do this inside the main loop because we also access the jth entry
	for(int i = 0; i < totalParticles; i++)
            for(j = 0; j < DIM; j++)
                  acceleration[i][j] = 0.0;

	for (int i = 0; i < totalParticles; i++)
	{
		for (j = 0; j < i; j++)
		{
                  pythagorean = 0;
                  for(k = 0; k < DIM; k++)
                  {
                        vectors[k] = determineVectorFlat(position[i][k], position[j][k]);
                        //vectors[k] = determineVectorPeriodic(position[i][k], position[j][k], boundaries[k]);
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

                        for(int k = 0; k < DIM; k++)
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
	for (int i = particlesType1; i < totalParticles; ++i)
            for(j = 0; j < DIM; j++)
                  acceleration[i][j] *= 0.5;
}

/*
void calcAcceleration(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles, 
	double particlesType1, double& potentialEnergy, double* boundaries, Tree &tree)
{
      //Zero out current acceleration
	for(int i = 0; i < totalParticles; i++)
	{
		acceleration[i][0] = 0.0;
		acceleration[i][1] = 0.0;
		acceleration[i][2] = 0.0;
            tree.calcAcc(0, i, position, acceleration[i])
	}
}
*/

void performEulerOperation(int totalParticles, double (*position)[DIM],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[DIM], double (*acceleration)[DIM],
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

                  applySolidBoundary(position[i][j], oldPosition[i][j], boundaries[j]);
                  //applyPeriodicBoundary(position[i][j], oldPosition[i][j], boundaries[j]);
			dotProd += velocity[i][j] * velocity[i][j];
		}

		kineticEnergy += (i < particlesType1) ? 0.5 * dotProd : dotProd;
	}
}

void performVerletOperation(int totalParticles, double (*position)[DIM],
	int particlesType1, double& potentialEnergy, double& kineticEnergy,
	double* boundaries, double (*oldPosition)[DIM], double (*acceleration)[DIM],
	double (*velocity)[DIM], double timestep, double halfInvTimestep)
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
        	for(j = 0; j < DIM; ++j) //Loop over all directions
		{
			currentDisplacement = position[i][j] - oldPosition[i][j];
			futureDisplacement = currentDisplacement + (dtsq * acceleration[i][j]);
			currentPosition = position[i][j];
			position[i][j] = position[i][j] + futureDisplacement;
			velocity[i][j] = (position[i][j] - oldPosition[i][j]) * (halfInvTimestep);
			oldPosition[i][j] = currentPosition;

                  applySolidBoundary(position[i][j], oldPosition[i][j], boundaries[j]);
                  //applyPeriodicBoundary(position[i][j], oldPosition[i][j], boundaries[j]);
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
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[DIM], const int totalParticles)
{	
      std::string str;
      //Fill out remaining zeros
      for(int i = 0; i < 3 - DIM; ++i)
            str.append("0 ");
      str.append("0 0 0 0 0 0 0 0 0\n");
	positionFile << "* " << currentTime << "\n";
	for (int i = 0; i < totalParticles; ++i)
      {
            for(int j = 0; j < DIM; ++j)
                  positionFile << position[i][j] << " ";
            positionFile << str;
      }
	//positionFile << position[i][0] << " " << position[i][1] << " " << position[i][2] << str;
}
