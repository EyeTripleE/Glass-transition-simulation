/*
*This program simulates the movement of inert particles using Lennard-Jones potential.
*Force is calculated from the potential.  Position and velocity are updated using the Verlet
* Algorithm (with Euler algorithm initialization).
*/

//NOTE: If running with OpenMP, make sure OMP_THREAD_LIMIT is set (e.g. export OMP_THREAD_LIMIT=5).
//This should probably be equal to OMP_NUM_THREADS.

//TODO: Add softening if ever use this again, figure out how to Speed up BH force
//calculation. Use 'restrict' keywork

//README!!!!!!!!!!!!! IMPORTANT: Set these environment variables before running
/*
export OMP_NUM_THREADS=<value>
export OMP_THREAD_LIMIT=<equal to OMP_NUM_THREADS>
export OMP_NESTED=TRUE

or on Intel Xeon Phi
export KMP_HOT_TEAMS_MODE=1
export KMP_HOT_TEAMS_MAX_LEVEL=2
export KMP_PLACE_THREADS=1T
export KMP_AFFINITY=compact,granularity=fine
*/

#include <cmath>
#include <iostream>
#include <ctime>
#include <fstream>
#include <random>
#include <vector>
#include <cstring>
#include "mpi.h"
#include <omp.h>
#include <chrono>
#include "assert.h"
#include <algorithm>

#define DIM 3

#define BARNES_HUT //Barnes Hut, tiles, or strips?
//#define TILES
#define CUTOFF //If strips, use cutoff distance or not?
#define OUTPUT //print output?

//===================================BEGIN FUNCTION HEADERS===================================
//Determines shortest vector from particle 1 to particle 2 (including across boundary) in one direction
inline double determineVectorPeriodic(const double p1Pos, const double p2Pos, const double size)
{
    double ds = p1Pos - p2Pos;
    return (fabs(ds) > 0.5 * size) ? (ds - copysign(size, ds)) : ds;
}

inline double determineVectorFlat(const double p1Pos, const double p2Pos)
{
    return p1Pos - p2Pos;
}

//===================================END FUNCTION HEADERS=====================================

//===================================BEGIN TREE STUFF=========================================

struct Node {
    float mass; // < 0 = LEAF, 0 = EMPTY_LEAF, > 0 = BRANCH
    double com[DIM];
    //double myBoundaries[6];
    //Run out of memory if save.
    //double boundaries[8][6];
    //std::vector<int> indicesArray[8];
};

//TODO Fix so works with 2D or 3D simulations
class Tree {

public:

    std::vector<Node> nodesArray;

    Tree()
    {

    }

    Tree(std::vector<Node> &vec)
    {
        nodesArray = vec;
    }

    void determineEntryPoints(int nodeIndex, double boundaries[6], int rank, int size,
                              std::vector<int> &entryNodes, double entryBoundaries[8][6])
    {

        double halves[3] = {0.5*(boundaries[0] + boundaries[1]),
                            0.5*(boundaries[2] + boundaries[3]),
                            0.5*(boundaries[4] + boundaries[5])
                           };

        double childBoundaries[8][6] = {{boundaries[0],halves[0],boundaries[2],halves[1],boundaries[4],halves[2]},//bottom, front, left
            {boundaries[0],halves[0],boundaries[2],halves[1],halves[2],boundaries[5]},
            {boundaries[0],halves[0],halves[1],boundaries[3],boundaries[4],halves[2]},
            {boundaries[0],halves[0],halves[1],boundaries[3],halves[2],boundaries[5]},
            {halves[0],boundaries[1],boundaries[2],halves[1],boundaries[4],halves[2]},
            {halves[0],boundaries[1],boundaries[2],halves[1],halves[2],boundaries[5]},
            {halves[0],boundaries[1],halves[1],boundaries[3],boundaries[4],halves[2]},
            {halves[0],boundaries[1],halves[1],boundaries[3],halves[2],boundaries[5]}
        };

        //Old code for calculating child boundaries
        /*for(int octant = 0; octant < 8; octant++)

        //Figure out the boundaries
        int index;
        int index2;
        int andVal = 7;
        int ltVal = 4; //less than value
        int i;

        {
              for(i = 0; i < 3; i++)
              {
                    index = 2*i;
                    index2 = index + 1;

                    if((octant & andVal) < ltVal)
                    {
                          entryBoundaries[octant][index] = boundaries[index];
                          entryBoundaries[octant][index2] = halves[i];
                    }
                    else
                    {
                          entryBoundaries[octant][index] = halves[i];
                          entryBoundaries[octant][index2] = boundaries[index2];
                    }

                    andVal >>= 1;
                    ltVal >>= 1;
              }
        }*/

        //Figure out the branches

        int level = floor( log(nodeIndex + 1) / log(8.0) ); //Get level of current node index
        int levelNodes = 1;
        for(int i = 0; i < level; i++) levelNodes *= 8;        
        int nextLevelNodes = levelNodes * 8;
        int index = (rank % nextLevelNodes) / levelNodes; //Convert rank to branch number

        if(size >= nextLevelNodes) //Too many processors, divide them to work on children
        {
            //Take next branch down the tree
            determineEntryPoints(8*nodeIndex + (index + 1), childBoundaries[index],
                                 rank, size, entryNodes, entryBoundaries);
        }
        else //No longer oversaturated
        {
            //This is only done once so it doesn't matter that unneeded information is copied
            memcpy(entryBoundaries, childBoundaries, 8*6*sizeof(double));
            //We have found the split point, recursion terminates     

            
            //How many processes are on my branch?    
            int pob = size / levelNodes; //Processes on branch
            if(index < size % levelNodes) pob++; //Distribute the remainder          
            
            //What is my rank relative to other processes on branch 
            int brank = rank / levelNodes;//branch rank

            int bpp = 8 / pob; //Branch per process
            int r = 8 % pob;

            int start;
            if(brank < r)
            {
               bpp++;
               start = brank*bpp;
            }
            else
            {
               start = r*(bpp + 1) + (brank - r)*bpp;
            }                    
            
            int entryIndex;
            for(int i = start + bpp - 1; i >= start; i--)
            {
                printf("%d\n", i);       
                entryIndex = 8*nodeIndex + (i + 1);
                entryNodes.push_back(entryIndex);            
            }           
            
            /*
            int stride = size / levelNodes;
            //Dole out the remainder
            if(index < size % levelNodes) stride += 1;
            int entryIndex;
            
            for(int i = 7 - index; i >= 0; i -= stride)
            {
                entryIndex = 8*nodeIndex + (i + 1);
                entryNodes.push_back(entryIndex);
            }
            */          
        }
    }

    void buildTree(int nodeIndex, double (*position)[DIM], std::vector<int> &partIndices, double boundaries[6]/*, std::vector<int> &keys*/)
    {
        //if(rank == 0) printf("here %d\n", nodeIndex);
        //printf("here %d\n", nodeIndex);
        //printf("Node index: %d, Number of particles: %ld, Boundaries: %g %g %g %g %g %g\n", nodeIndex, partIndices.size(),
        //      boundaries[0], boundaries[1], boundaries[2],boundaries[3],boundaries[4],boundaries[5]);

        //If there aren't any writers, increment resizers
        //#pragma omp critical
        //{
        //   while(numWriters != 0 && numResizers != 0){printf("here1\n");};//Wait until the number of writers and resizers is zero
        //   numResizers++;  //Increment resizers.
        //}
        //If another process resizes while another thread accesses, program will crash

        #pragma omp critical
        if(nodeIndex >= nodesArray.size()) nodesArray.resize(nodeIndex + 1);

        //numResizers--; //Decrement resizers. Shouldn't need to be atomic because only one thread can modify at a time

        if(partIndices.size() == 0)
        {
            /*
            #pragma omp critical
            {               
               while(numResizers != 0){printf("here2\n");}; //Wait until the number of resizers is zero
               #pragma omp atomic
               numWriters++;  //Increment writers. These two operations must happen atomically
            }
            */

            #pragma omp critical
            nodesArray[nodeIndex].mass = 0.0f;

            /*
            #pragma omp atomic
            numWriters--; //Decrement writers
            */
        }
        else if(partIndices.size() == 1)
        {
            /*
            #pragma omp critical //This basically makes it serial, might as well just 
            //do a critical section around the writes if doing this.
            {               
               while(numResizers != 0){
               printf("here3\n");}; //Wait until the number of resizers is zero
               #pragma omp atomic
               numWriters++;  //Increment writers. These two operations must happen atomically
            }
            */
            #pragma omp critical
            {
            nodesArray[nodeIndex].mass = -1.0f; //Assume equivalent mass, otherwise will need to modify
            memcpy(&(nodesArray[nodeIndex].com[0]), &position[partIndices[0]][0], DIM*sizeof(double));
            //keys.push_back(partIndices[0]);
            }
            /*
            #pragma omp atomic
            numWriters--; //Decrement writers
            */
        }
        else
        {
            //Assuming mass of one per particle, average location
            //Calculate the center of mass
            double com0 = 0.0;
            double com1 = 0.0;
            double com2 = 0.0;

            //With OpenMP 4.0 an array (com[3]) could be reduced. Until then...
            #pragma omp parallel for reduction(+:com0, com1, com2)
            for(int i = 0; i < partIndices.size(); i++)
            {
                com0 += position[partIndices[i]][0];
                com1 += position[partIndices[i]][1];
                com2 += position[partIndices[i]][DIM - 1]; //Same as com1 if DIM == 2, messes up vectorization for two dimensions
            }
            double invSize = 1.0/partIndices.size();

            /*
            #pragma omp critical
            {               
               while(numResizers != 0){
               printf("here4\n");
               }; //Wait until the number of resizers is zero
               #pragma omp atomic
               numWriters++;  //Increment writers. These two operations must happen atomically
            }
            */

            #pragma omp critical
            {
            nodesArray[nodeIndex].com[0] = com0*invSize;
            nodesArray[nodeIndex].com[1] = com1*invSize;
            nodesArray[nodeIndex].com[DIM - 1] = com2*invSize;
            nodesArray[nodeIndex].mass = partIndices.size(); //Assuming mass of one, sum mass
            }

            /*
            #pragma omp atomic
            numWriters--; //Decrement writers
            */

            /*
            //Old code for calculating COM
            //Zero out center of mass
            memset(nodesArray[nodeIndex].com, 0, DIM*sizeof(double));
            
            int j;
            for(int i = 0; i < partIndices.size(); i++)
            {
                  index = partIndices[i];
                  for(j = 0; j < DIM; j++)
                  {
                        nodesArray[nodeIndex].com[j] += position[index][j];
                  }
            }            

            double invSize = 1.0/partIndices.size();
            for(int i = 0; i < DIM; i++)
            {
                nodesArray[nodeIndex].com[i] *= invSize;
            }
            */

            //Create new boundaries
            double halves[3] = {0.5*(boundaries[0] + boundaries[1]),
                                0.5*(boundaries[2] + boundaries[3]),
                                0.5*(boundaries[4] + boundaries[5])
                               };

            double childBoundaries[8][6] = {{boundaries[0],halves[0],boundaries[2],halves[1],boundaries[4],halves[2]},//bottom, front, left
                {boundaries[0],halves[0],boundaries[2],halves[1],halves[2],boundaries[5]},
                {boundaries[0],halves[0],halves[1],boundaries[3],boundaries[4],halves[2]},
                {boundaries[0],halves[0],halves[1],boundaries[3],halves[2],boundaries[5]},
                {halves[0],boundaries[1],boundaries[2],halves[1],boundaries[4],halves[2]},
                {halves[0],boundaries[1],boundaries[2],halves[1],halves[2],boundaries[5]},
                {halves[0],boundaries[1],halves[1],boundaries[3],boundaries[4],halves[2]},
                {halves[0],boundaries[1],halves[1],boundaries[3],halves[2],boundaries[5]}
            };

            std::vector<int> indicesArray[8];

            //Subdivide the indices based on boundaries
            //push_back means OpenMP cannot perform well
            int childOctant;
            for(int i = 0; i < partIndices.size(); i++)
            {
                //Have to zero out each iteration if using if statement
                childOctant = position[partIndices[i]][0] > halves[0] ? 4 : 0;
                if(position[partIndices[i]][1] > halves[1]) childOctant |= 2;
                if(position[partIndices[i]][2] > halves[2]) childOctant |= 1;
                indicesArray[childOctant].push_back(partIndices[i]);
            }

            #pragma omp parallel for
            for(int i = 7; i >= 0; i--) //Reduce the number of resizes by going right to left
            //for(int i = 0; i < 8; i++)            
            {
                buildTree(8*nodeIndex + (i + 1), position, indicesArray[i], childBoundaries[i]/*, keys*/);
            }
        }
    }

    void calcAcc(int nodeIndex, std::vector<int> &indices, double (*position)[DIM], double (*acceleration)[DIM], double &potentialEnergy)
    {
        //Ensure we are not at bottom of tree and we still have particles to work with
        if(indices.size() > 0 && nodesArray[nodeIndex].mass != 0.0f)
        {
            //Create vector for the next recursion level
            std::vector<int> indicesForChildren;

            //Initialize variables for use in loop
            double vectors[DIM];
            double sigma, sigmaPow6, sigmaPow12;
            double invPy, invPyPow3, invPyPow4, invPyPow6;
            double forceCoeff, pythagorean, force;
            double cutoff_sq = 5.0;//1.0*1.0;//0.5*0.5;
            float mass_sq = nodesArray[nodeIndex].mass*nodesArray[nodeIndex].mass; 
            double m_over_c = mass_sq / cutoff_sq;
            int j, k, index;

            //Look at all my particles
            for(int i = 0; i < indices.size(); i++) //This is difficult to vectorize unless particles are arranged in spatial order
            {
               index = indices[i];
               pythagorean = 0.0;
               for(j = 0; j < DIM; j++)
               {
                   //This random access is likely costly, but sorting particles doesn't seem to reduce gap sizes 
                   vectors[j] = determineVectorFlat(position[index][j], nodesArray[nodeIndex].com[j]);
                   pythagorean += vectors[j]*vectors[j];
               }             

               //With higher s values may need to decrease timestep. If pythagorean == 0.0 then particle is looking at itself
               if(m_over_c < pythagorean || (nodesArray[nodeIndex].mass < 0.0f && pythagorean != 0.0)) //Calculate force
               {
                   //Guessing this should be an average of sigma values of all involved particles. Since only using type 1 just set to one here.
                   sigma = 1.0;//(i < particlesType1 && j < particlesType1) ? 1.0 : ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);
                   sigmaPow6 = sigma*sigma;//Use sigmaPow6 as to hold the square temporarily
                   sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
                   sigmaPow12 = sigmaPow6*sigmaPow6;
                   invPy = 1.0 / pythagorean;
                   invPyPow3 = invPy*invPy*invPy;
                   invPyPow4 = invPyPow3*invPy;
                   invPyPow6 = invPyPow3*invPyPow3;
                   potentialEnergy += 2.0 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
                   forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);

		   for(k = 0; k < DIM; k++)
                   {
                       force = vectors[k]*forceCoeff;
                       #pragma omp atomic update //Sibling nodes could try to update
                       acceleration[index][k] += force; //Another random access
                   }
               }
               else if(pythagorean != 0.0) //Add to vectors for child nodes
               {
                  indicesForChildren.push_back(index);
               }
            }

            //Recurse to children
            double childPE = 0.0;
            #pragma omp parallel for reduction(+:childPE)
            for(int i = 0; i < 8; i++)
               calcAcc(8*nodeIndex + (i + 1), indicesForChildren, position, acceleration, childPE);

            potentialEnergy += childPE;            
        }
    }

};

//============================END TREE STUFF============================================

//============================MORE FUNCTION DECLARATIONS======================================

//Calculates force on each particle based on the present position of all particles
//Cutoff distance of 4.0 units.
//Then calculates the acceleration of each particle based on the force
//currently applied to it.
void calcAccelerationStrips(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles,
                            int myStart, int myEnd, int particlesType1, double &potentialEnergy, double boundaries[6],
                            MPI_Request &request);

//Tiled version
void calcAccelerationTiles(double (*acceleration)[DIM], double (*position)[DIM],
                           int myStart, int myEnd, int rowStart, int rowEnd, int columnStart, int columnEnd,
                           int targetRank, int coord[2], int particlesType1, double &potentialEnergy, double boundaries[6],
                           MPI_Request &request);

//Barnes-Hut version
void calcAccelerationBH(double (*acceleration)[DIM], double (*position)[DIM], double (*oldPosition)[DIM],
                        int myStart, int myEnd, int rank, int size, int particlesType1, double &potentialEnergy,
                        Tree &tree, std::vector<int> indices[8], std::vector<int> &allIndices, std::vector<int> &entryNodes,
                        double entryBoundaries[8][6], int entryOctants[8], MPI_Request &request);

//Function that performs the Euler algorithm on all particles in the set
void performEulerOperation(int myStart, int myEnd, int accDisp,
                           double (*position)[DIM], double (*oldPosition)[DIM], double (*acceleration)[DIM], double (*velocity)[DIM],
                           double boundaries[6], double timestep, MPI_Request &request);

//Function that performs the Verlet algorithm on all particles in the set
void performVerletOperation(int myStart, int myEnd, int accDisp, double (*position)[DIM], double (*oldPositionHolder)[DIM],
                            double (*oldPosition)[DIM], double (*acceleration)[DIM], double boundaries[6], double dtsq, MPI_Request &request);

//Calculates the kinetic energy at the previous time step
double calcKineticEnergy(int myStart, int myEnd, double (*position)[DIM], double (*oldPosition)[DIM],
                         int particlesType1, double invTimestep);

//Translates particles that have exited the simulation area back into the
//Simulation region, works in one direction so do for every direction
void applyPeriodicBoundary(double &position, double &oldPosition, double boundary[2]);

//Bounce particle off wall
void applySolidBoundary(double &position, double &oldPosition, double boundary[2]);

//Outputs to position file
void outputPosition(std::ofstream &positionFile, const double currentTime, double (*position)[DIM], const int totalParticles);

//Initializes
int initializeParticles(double (*position)[DIM], double (*velocity)[DIM],
                        int totalParticles, double boundaries[6]);

//Deletes arrays and closes files
void cleanup(double (*position)[DIM], double (*velocity)[DIM], double (*acceleration)[DIM],
             double (*oldPosition)[DIM]);

//=======================================================================================

MPI_Comm COMM_COLUMN;
MPI_Comm COMM_ROW;

//The main function to execute the simulation
int main(int argc, char* argv[])
{
#if defined(TILES) && defined(BARNES_HUT)
    printf("Error: Incompatible precompiler directives\n");
    exit(1);
#endif

    //Setup MPI
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Request request;

#if defined(TILES)
    if(!(size == 1 || size == 2 || size == 4 || size == 16 || size == 64))
    {
        printf("Error: Number of processes is not a (small) factor of 1024\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif

    //START SETUP
    //clock_t tstart = clock();
    /*CONSTANTS FOR REFERENCE THESE ARE HARDCODED*/
    //constants 0.5(sigma1 + sigma2)
    //double sigma1to1 = 1; //in units of sigma
    //double sigma2to2 = 1.4; //in units of sigma
    //double sigma1to2 = 1.2; // in units of sigma
    //double epsilon = 1;// in units of epsilon
    //double massParticle1 = 1; //in units of massParticle1
    //double massParticle2 = 2; //in units of massParticle1

    //Set variables for number of particles
    int numParticlesType1 = 1024;
    if (argc > 1)
        numParticlesType1 = atoi(argv[1]);
    int numParticlesType2 = 0;
    int totalParticles = numParticlesType1 + numParticlesType2;

    //Set MPI variables
    int localParticles = totalParticles/size;
    int myStart = rank*localParticles;
    int myEnd = (rank + 1)*localParticles;

#if defined(TILES)

    int reorder = 0;
    int coord[2];
    int dim[2] = {0, 0};
    int period[2] = {0, 0};

    MPI_Dims_create(size, 2, dim);
    assert(dim[0] == dim[1]);

    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm_cart);
    MPI_Cart_coords(comm_cart, rank, 2, coord);

    int groupParticles = totalParticles / dim[0];
    int rowStart = coord[0]*groupParticles;
    int rowEnd = (coord[0] + 1)*groupParticles;
    int columnStart = coord[1]*groupParticles;
    int columnEnd = (coord[1] + 1)*groupParticles;

    int target2D[2] = {coord[1], coord[0]};
    int targetRank;
    MPI_Cart_rank(comm_cart, target2D, &targetRank);

    // dim[0] is x and dim[1] is y
    MPI_Comm_split(MPI_COMM_WORLD, coord[0], rank, &COMM_ROW);
    MPI_Comm_split(MPI_COMM_WORLD, coord[1], rank, &COMM_COLUMN);

#endif

    //Set time related variables
    double timestep = 0.00005; //Can be arbitrarily small
    double maxTime = 1000*timestep; //Can be arbitrarily long or short
    double invTimestep = 1.0/timestep; //needed for KE
    double dtsq = timestep*timestep;
    double currentTime = 0;

    //Set energy variables
    double energies[2], totalEnergy; //{kineticEnergy, potentialEnergy}

    //Set boundary variables
    double boundaries[6] = {0.0, 50.0, 0.0, 50.0, 0.0, 50.0}; //{x_min, x_max, y_min, y_max, z_min, z_max}

    //Particle information arrays
    double (*position)[DIM] = new double[totalParticles][DIM];
    double (*velocity)[DIM] = new double[totalParticles][DIM];
    double (*oldPosition)[DIM] = new double[totalParticles][DIM];
    double (*acceleration)[DIM] = new double[totalParticles][DIM];
    double (*oldPositionHolder)[DIM] = velocity; //Reuse the location for the oldPositionHolder

    int errorStatus = 0;

    if(rank == 0)
    {
        errorStatus = initializeParticles(position, velocity, totalParticles, boundaries);
    }

    MPI_Bcast(&errorStatus, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(errorStatus != 0)
    {
        cleanup(position, velocity, acceleration, oldPosition);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //-----------------------------------------------------------------------//
    // Distribute Data           						             //
    MPI_Bcast(position, DIM*totalParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);	 //
    MPI_Bcast(velocity, DIM*totalParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);	 //
    //-----------------------------------------------------------------------//

    //Create data files
    std::ofstream positionFile, energyFile;
    if(rank == 0)
    {
        positionFile.open("position.txt");
        energyFile.open("energy.txt");
    }

#if defined(TILES)
    int accDisp = rowStart - myStart;
#elif defined(BARNES_HUT)
    int accDisp = -myStart; //Start at zero
    Tree tree;
    double entryBoundaries[8][6];
    std::vector<int> entryNodes;
    tree.determineEntryPoints(0, boundaries, rank, size, entryNodes, entryBoundaries);
    std::vector<int> indices[8];
    int entryOctants[8];
    for(int i = 0; i < entryNodes.size(); i++)
        entryOctants[i] = entryNodes[i] - ((((entryNodes[i] - 1) / 8 ) * 8) + 1);
    std::vector<int> allIndices;
    for(int i = 0; i < totalParticles; i++)
        allIndices.push_back(i);
#else
    int accDisp = 0;
#endif

    //END SETUP

    //START SIMULATION
    auto start_time = std::chrono::high_resolution_clock::now(); //Start the timer

#if defined(BARNES_HUT)
    calcAccelerationBH(acceleration, position, oldPosition, myStart, myEnd, rank, size, numParticlesType1, energies[1],
                       tree, indices, allIndices, entryNodes, entryBoundaries, entryOctants, request);
#elif defined(TILES)
    calcAccelerationTiles(acceleration, position, myStart, myEnd, rowStart, rowEnd, columnStart,
                          columnEnd, targetRank, coord, numParticlesType1, energies[1], boundaries, request);
#else
    calcAccelerationStrips(acceleration, position, totalParticles, myStart, myEnd,
                           numParticlesType1, energies[1], boundaries, request);
#endif

    //Perform initial Euler operation to set things in motion, actually a second order Taylor,
    //not first order Euler.
    performEulerOperation(myStart, myEnd, accDisp, position, oldPosition, acceleration, velocity,
                          boundaries, timestep, request);

    //Main loop - performing Verlet operations for the remainder of the simulation
    unsigned count = 0;
    for (currentTime = 2*timestep; currentTime < maxTime; currentTime += timestep)
    {
#if defined(BARNES_HUT)
        calcAccelerationBH(acceleration, position, oldPosition, myStart, myEnd, rank, size,
                           numParticlesType1, energies[1], tree, indices, allIndices, entryNodes,
                           entryBoundaries, entryOctants, request);
#elif defined(TILES)
        calcAccelerationTiles(acceleration, position, myStart, myEnd, rowStart, rowEnd,
                              columnStart, columnEnd, targetRank, coord, numParticlesType1, energies[1],
                              boundaries, request);
#else
        calcAccelerationStrips(acceleration, position, totalParticles, myStart, myEnd,
                               numParticlesType1, energies[1], boundaries, request);
#endif

        performVerletOperation(myStart, myEnd, accDisp, position, oldPositionHolder,
                               oldPosition, acceleration, boundaries, dtsq, request);

#ifdef OUTPUT
        count++;  //Can set print interval arbitrarily
        if(count >= 10)
        {
            count = 0;
            energies[0] = calcKineticEnergy(myStart, myEnd, position, oldPosition, numParticlesType1, invTimestep);
            if(rank == 0)
            {
                MPI_Reduce(MPI_IN_PLACE, energies, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

                //Ideally would spin off thread to handle, since we're benchmarking without OUTPUT
                //this isn't necessary.
                totalEnergy = energies[0] + energies[1];
                printf("%g %g %g %g\n", currentTime, energies[0], energies[1], totalEnergy);
                outputPosition(positionFile, currentTime, position, totalParticles);
                energyFile << currentTime << " " << totalEnergy << " " << energies[0]
                           << " " << energies[1] << "\n";
            }
            else
            {
                MPI_Reduce(energies, nullptr, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
        }
#endif
    }

    cleanup(position, velocity, acceleration, oldPosition);

    if(rank == 0)
    {
        //printf("Total time (s): %g\n", (double)(clock() - tstart) / CLOCKS_PER_SEC);
        std::chrono::high_resolution_clock::duration diff = std::chrono::high_resolution_clock::now()-start_time;
        std::cout<<"Ranks = " << size << " Time taken = "<<diff.count()<<std::endl;
        positionFile.close();
        energyFile.close();
    }


    MPI_Finalize();

    return 0;
}
//End of main function

void calcAccelerationStrips(double (*acceleration)[DIM], double (*position)[DIM], double totalParticles,
                            int myStart, int myEnd, int particlesType1, double &potentialEnergy, double boundaries[6], MPI_Request &request)
{
    int localSize = DIM*(myEnd - myStart);
    MPI_Iallgather(MPI_IN_PLACE, localSize, MPI_DOUBLE, position,
                   localSize, MPI_DOUBLE, MPI_COMM_WORLD, &request);

    double sigma, sigmaPow6, sigmaPow12;
    double pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6;
    double vectors[DIM], forceCoeff;
    int j, k;

    memset(&acceleration[myStart], 0, localSize*sizeof(double));
    double pe = 0;

    MPI_Wait(&request, MPI_STATUS_IGNORE);

    #pragma omp parallel for reduction(+:pe) private(j, k, sigma, sigmaPow6, sigmaPow12, \
    pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6, vectors, forceCoeff)
    for (int i = myStart; i < myEnd; i++)
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

#ifdef CUTOFF
            if (pythagorean < 16.0)
#endif
            {
                //Force derived from Lennard-Jones potential
                sigma = 1.0;//(i < particlesType1 && j < particlesType1) ? 1.0 : ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);
                sigmaPow6 = sigma*sigma;//Use sigmaPow6 to hold the square temporarily
                sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
                sigmaPow12 = sigmaPow6*sigmaPow6;
                invPy = 1.0 / pythagorean;
                invPyPow3 = invPy*invPy*invPy;
                invPyPow4 = invPyPow3*invPy;
                invPyPow6 = invPyPow3*invPyPow3;
                pe += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
                forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);

                for(k = 0; k < DIM; k++)
                {
                    acceleration[i][k] += vectors[k]*forceCoeff;
                }
            }
        }
        for (j = i + 1; j < totalParticles; j++)
        {
            pythagorean = 0;
            for(k = 0; k < DIM; k++)
            {
                vectors[k] = determineVectorFlat(position[i][k], position[j][k]);
                //vectors[k] = determineVectorPeriodic(position[i][k], position[j][k], boundaries[2*k + 1]);
                pythagorean += vectors[k]*vectors[k];
            }

#ifdef CUTOFF
            if (pythagorean < 16.0)
#endif
            {
                //Force derived from Lennard-Jones potential
                sigma = 1.0;//(i < particlesType1 && j < particlesType1) ? 1.0 : ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);
                sigmaPow6 = sigma*sigma;//Use sigmaPow6 to hold the square temporarily
                sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
                sigmaPow12 = sigmaPow6*sigmaPow6;
                invPy = 1.0 / pythagorean;
                invPyPow3 = invPy*invPy*invPy;
                invPyPow4 = invPyPow3*invPy;
                invPyPow6 = invPyPow3*invPyPow3;
                pe += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
                forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);

                for(k = 0; k < DIM; k++)
                {
                    acceleration[i][k] += vectors[k]*forceCoeff;
                }
            }
        }
    }

    potentialEnergy = pe;

    //Because type two particles are twice as heavy
    //for (i = particlesType1; i < totalParticles; ++i)
    //      for(j = 0; j < DIM; j++)
    //            acceleration[i][j] *= 0.5;
}


void calcAccelerationTiles(double (*acceleration)[DIM], double (*position)[DIM],
                           int myStart, int myEnd, int rowStart, int rowEnd, int columnStart, int columnEnd,
                           int targetRank, int coord[2], int particlesType1, double &potentialEnergy, double boundaries[6], MPI_Request &request)
{
    int localSize = DIM*(myEnd - myStart);
    int groupSize = DIM*(rowEnd - rowStart);
    MPI_Request temp;

    //Gather along the rows
    MPI_Allgather(MPI_IN_PLACE, localSize, MPI_DOUBLE, &position[rowStart],
                  localSize, MPI_DOUBLE, COMM_ROW);

    //Point to point communication, it's a transpose!
    if(coord[0] != coord[1])
    {
        MPI_Isend(&position[rowStart], groupSize, MPI_DOUBLE,
                  targetRank, 0, MPI_COMM_WORLD, &temp);

        MPI_Irecv(&position[columnStart], groupSize, MPI_DOUBLE,
                  MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
    }

    double sigma, sigmaPow6, sigmaPow12;
    double pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6;
    double vectors[DIM], forceCoeff;
    int j, k;

    memset(&acceleration[rowStart], 0, groupSize*sizeof(double));
    double pe = 0.0;


    if(coord[0] != coord[1])
    {
        //Don't need to wait for the send to finish
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    //Calculate force
    if(coord[0] == coord[1])
    {
        #pragma omp parallel for reduction(+:pe) private(j, k, sigma, sigmaPow6, sigmaPow12, \
        pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6, vectors, forceCoeff)
        for (int i = rowStart; i < rowEnd; i++)
    {
        for (j = columnStart; j < i; j++)
            {
                pythagorean = 0;
                for(k = 0; k < DIM; k++)
                {
                    vectors[k] = determineVectorFlat(position[i][k], position[j][k]);
                    //vectors[k] = determineVectorPeriodic(position[i][k], position[j][k], boundaries[2*k + 1]);
                    pythagorean += vectors[k]*vectors[k];
                }

#ifdef CUTOFF
                if (pythagorean < 16.0)
#endif
                {
                    //Force derived from Lennard-Jones potential
                    sigma = 1.0;//(i < particlesType1 && j < particlesType1) ? 1.0 : ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);
                    sigmaPow6 = sigma*sigma;//Use sigmaPow6 to hold the square temporarily
                    sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
                    sigmaPow12 = sigmaPow6*sigmaPow6;
                    invPy = 1.0 / pythagorean;
                    invPyPow3 = invPy*invPy*invPy;
                    invPyPow4 = invPyPow3*invPy;
                    invPyPow6 = invPyPow3*invPyPow3;
                    pe += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
                    forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);

                    for(k = 0; k < DIM; k++)
                    {
                        acceleration[i][k] += vectors[k]*forceCoeff;
                    }
                }
            }
            for (j = i + 1; j < columnEnd; j++)
            {
                pythagorean = 0;
                for(k = 0; k < DIM; k++)
                {
                    vectors[k] = determineVectorFlat(position[i][k], position[j][k]);
                    //vectors[k] = determineVectorPeriodic(position[i][k], position[j][k], boundaries[2*k + 1]);
                    pythagorean += vectors[k]*vectors[k];
                }

#ifdef CUTOFF
                if (pythagorean < 16.0)
#endif
                {
                    //Force derived from Lennard-Jones potential
                    sigma = 1.0;//(i < particlesType1 && j < particlesType1) ? 1.0 :
                    //((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);
                    sigmaPow6 = sigma*sigma;//Use sigmaPow6 to hold the square temporarily
                    sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
                    sigmaPow12 = sigmaPow6*sigmaPow6;
                    invPy = 1.0 / pythagorean;
                    invPyPow3 = invPy*invPy*invPy;
                    invPyPow4 = invPyPow3*invPy;
                    invPyPow6 = invPyPow3*invPyPow3;
                    pe += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
                    forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);

                    for(k = 0; k < DIM; k++)
                    {
                        acceleration[i][k] += vectors[k]*forceCoeff;
                    }
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for reduction(+:pe) private(j, k, sigma, sigmaPow6, sigmaPow12, \
        pythagorean, invPy, invPyPow3, invPyPow4, invPyPow6, vectors, forceCoeff)
        for (int i = rowStart; i < rowEnd; i++)
    {
        for (j = columnStart; j < columnEnd; j++)
            {
                pythagorean = 0;
                for(k = 0; k < DIM; k++)
                {
                    vectors[k] = determineVectorFlat(position[i][k], position[j][k]);
                    //vectors[k] = determineVectorPeriodic(position[i][k], position[j][k], boundaries[2*k + 1]);
                    pythagorean += vectors[k]*vectors[k];
                }

#ifdef CUTOFF
                if (pythagorean < 16.0)
#endif
                {
                    //Force derived from Lennard-Jones potential

                    sigma = 1.0;//(i < particlesType1 && j < particlesType1) ? 1.0 :
                    // ((i >= particlesType1 && j >= particlesType1) ? 1.4 : 1.2);
                    sigmaPow6 = sigma*sigma;//Use sigmaPow6 to hold the square temporarily
                    sigmaPow6 = sigmaPow6*sigmaPow6*sigmaPow6;
                    sigmaPow12 = sigmaPow6*sigmaPow6;
                    invPy = 1.0 / pythagorean;
                    invPyPow3 = invPy*invPy*invPy;
                    invPyPow4 = invPyPow3*invPy;
                    invPyPow6 = invPyPow3*invPyPow3;
                    pe += 2 * ((sigmaPow12 * invPyPow6) - (sigmaPow6 * invPyPow3));
                    forceCoeff = (sigmaPow6 * invPyPow4) * ((48.0 * sigmaPow6 * invPyPow3) - 24.0);

                    for(k = 0; k < DIM; k++)
                    {
                        acceleration[i][k] += vectors[k]*forceCoeff;
                    }
                }
            }
        }
    }
    //Reduce and scatter, send to other processors in same row
    MPI_Ireduce_scatter_block(MPI_IN_PLACE, &acceleration[rowStart],
                              localSize, MPI_DOUBLE, MPI_SUM, COMM_ROW, &request);

    potentialEnergy = pe;
}

void calcAccelerationBH(double (*acceleration)[DIM], double (*position)[DIM], double (*oldPosition)[DIM],
                        int myStart, int myEnd, int rank, int size, int particlesType1, double &potentialEnergy,
                        Tree &tree, std::vector<int> indices[8], std::vector<int> &allIndices, std::vector<int> &entryNodes,
                        double entryBoundaries[8][6], int entryOctants[8], MPI_Request &request)
{
    int localParticles = myEnd - myStart;
    int localSize = DIM*localParticles;

    //MPI_Request requests[2];
    //This is expensive, need to get rid of
    MPI_Iallgather(MPI_IN_PLACE, localSize, MPI_DOUBLE, position,
                   localSize, MPI_DOUBLE, MPI_COMM_WORLD, &request);
    //Only one process has the updated oldPosition information. 
    //MPI_Iallgather(MPI_IN_PLACE, localSize, MPI_DOUBLE, oldPosition,
    //               localSize, MPI_DOUBLE, MPI_COMM_WORLD, &requests[1]);

    int totalParticles = size*localParticles;
    int totalSize = localSize*size;

    //Set acceleration of particles to zero, this is legal
    //http://stackoverflow.com/questions/4629853/is-it-legal-to-use-memset-0-on-array-of-doubles
    memset(acceleration, 0, totalSize*sizeof(double));
    double pe = 0;
    int j, k;

    for(int i = 0; i < entryNodes.size(); i++)
    {
        indices[i].clear();
    }

    MPI_Wait(&request, MPI_STATUS_IGNORE);
    //MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

    //Maybe overlap the tree building and communication? How?
    //Not necessary. Tree building is not a significant cost.

    //Loop over particles and add it to the corresponding list of particles
    //if it is inside one of my octants
    #pragma omp parallel for private(j, k)
    for(int i = 0; i < totalParticles; i++)
    {
        for(j = 0; j < entryNodes.size(); j++)
        {
            for(k = 0; k < 3 && position[i][k] > entryBoundaries[entryOctants[j]][2*k] &&
                        position[i][k] <= entryBoundaries[entryOctants[j]][2*k + 1]; k++)
            {

            }

            if(k == 3)
            {
                #pragma omp critical
                indices[j].push_back(i);

                break;
            }
        }
    }

   /*
   std::vector<int> keys;
   keys.reserve(totalParticles);
   std::vector<int> recvCounts(size);
   std::vector<int> disp(size); disp[0] = 0;
   std::vector<int> globalKeys(totalParticles);
   */

   #pragma omp parallel for //Not thread safe due to overlapping resizes and writes
   for(int i = 0; i < entryNodes.size(); i++)
   {
      tree.buildTree(entryNodes[i], position, indices[i], entryBoundaries[entryOctants[i]]/*, keys*/);
   }

   //Sort particles so that the particle at the index of the key's value goes to that position
   //Basic version: Send counts to each process, then send keys, then sort entire list locally
   
   /*
   recvCounts[rank] = keys.size();
   MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, &recvCounts[0], 1, MPI_INT, MPI_COMM_WORLD);

   for(int i = 1; i < size; i++)
      disp[i] = disp[i - 1] + recvCounts[i - 1];

   MPI_Allgatherv(&keys[0], keys.size(), MPI_INT, &globalKeys[0], &recvCounts[0], &disp[0], MPI_INT, MPI_COMM_WORLD);

   //Why is it completely different each time? Shouldn't consecutive trees
   //be nearly identical?  
   */
   /*
   if(rank == 0)
   {
      for(int i = 0; i < totalParticles; i++)
      {
         printf("%d %d\n", globalKeys[i], allIndices[i]);
      }
      printf("\n\n");
   } 
   */
   

   //More advanced version - make size bins of width n/size and bin local particles. Send particles in a bin to
   //Appropriate processor and have each process sort all particles in a bin.
   //Then gather the results

   //Put particles data in correct location.
   
   /*
   double temp[DIM];
   int index;
   int s = DIM*sizeof(double);
   for(int i = 0; i < totalParticles - 1; i++)
   {
      index = globalKeys[i];
      while(i != index)
      {
         //Swap data, data at index will be in correct location
         memcpy(&temp[0], &position[index][0], s);
         memcpy(&position[index][0], &position[i][0], s);
         memcpy(&position[i][0], &temp[0], s);

         memcpy(&temp[0], &oldPosition[index][0], s);
         memcpy(&oldPosition[index][0], &oldPosition[i][0], s);
         memcpy(&oldPosition[i][0], &temp[0], s);
         
         //Update keys
         globalKeys[i] = globalKeys[index];
         globalKeys[index] = index;
         index = globalKeys[i];             
      }
   }   
   */

   /*
   if(rank == 0)
   {
      for(int i = 0; i < totalParticles; i++)
      {
         printf("%d %d\n", globalKeys[i], allIndices[i]);
      }
      printf("\n\n");
   } 
   */
   
   //Most time spent here.
   #pragma omp parallel for reduction(+:pe) //Can't do in above loop because tree resizes 
   for(int i = 0; i < entryNodes.size(); i++)
   {
      tree.calcAcc(entryNodes[i], allIndices, position, acceleration, pe);
   }

   //Is the problem that not all ranks have updated oldPosition?
   MPI_Ireduce_scatter_block(MPI_IN_PLACE, acceleration,
                           localSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request);

   potentialEnergy = pe;
}

//Function that performs the Euler algorithm on all particles in the set
void performEulerOperation(int myStart, int myEnd, int accDisp,
                           double (*position)[DIM], double (*oldPosition)[DIM], double (*acceleration)[DIM], double (*velocity)[DIM],
                           double boundaries[6], double timestep, MPI_Request &request)
{
    int j;
    double dtsq = timestep*timestep;

    memcpy(&oldPosition[myStart][0], &position[myStart][0], DIM*sizeof(double)*(myEnd - myStart));

#if defined(BARNES_HUT) || defined (TILES)
    MPI_Wait(&request, MPI_STATUS_IGNORE);
#endif

    #pragma omp parallel for private(j)
    for (int i = myStart; i < myEnd; i++)
    {
        for (j = 0; j < DIM; ++j)
        {
            position[i][j] += velocity[i][j] * timestep + 0.5*dtsq*acceleration[accDisp + i][j];
            applySolidBoundary(position[i][j], oldPosition[i][j], &boundaries[2*j]);
            //applyPeriodicBoundary(position[i][j], oldPosition[i][j], boundaries[2*j]);
        }
    }
}

void performVerletOperation(int myStart, int myEnd, int accDisp, double (*position)[DIM], double (*oldPositionHolder)[DIM],
                            double (*oldPosition)[DIM], double (*acceleration)[DIM], double boundaries[6], double dtsq, MPI_Request &request)
{
    int j;
    int numBytes = DIM*sizeof(double)*(myEnd - myStart);

    memcpy(&oldPositionHolder[myStart][0], &oldPosition[myStart][0], numBytes);
    memcpy(&oldPosition[myStart][0], &position[myStart][0], numBytes);

#if defined(BARNES_HUT) || defined (TILES)
    MPI_Wait(&request, MPI_STATUS_IGNORE);
#endif

    #pragma omp parallel for private(j)
    for (int i = myStart; i < myEnd; i++)
    {
        // Vector Verlet Method
        for(j = 0; j < DIM; ++j) //Loop over all directions
        {
            position[i][j] += (position[i][j] - oldPositionHolder[i][j]) + (dtsq * acceleration[accDisp + i][j]);
            applySolidBoundary(position[i][j], oldPosition[i][j], &boundaries[2*j]);
            //applyPeriodicBoundary(position[i][j], oldPosition[i][j], boundaries[2*j]);
        }
    }
}

double calcKineticEnergy(int myStart, int myEnd, double (*position)[DIM], double (*oldPosition)[DIM],
                         int particlesType1, double invTimestep)
{
    double velocity;
    double dotProd;
    double ke = 0;
    int j;
    #pragma omp parallel for reduction(+:ke) private(j, velocity, dotProd)
    for(int i = myStart; i < myEnd; i++)
    {
        dotProd = 0;
        for(j = 0; j < DIM; ++j)
        {
            velocity = (position[i][j] - oldPosition[i][j]) * invTimestep;
            dotProd += velocity*velocity;
        }
        ke += (i < particlesType1) ? 0.5 * dotProd : dotProd;
    }
    return ke;
}

//CHANGE THESE FUNCTIONS IF WANT TO SUPPORT NON-ZERO LOWER BOUND
void applyPeriodicBoundary(double &position, double &oldPosition, double boundary[2])
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

void applySolidBoundary(double &position, double &oldPosition, double boundary[2])
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

int initializeParticles(double (*position)[DIM], double (*velocity)[DIM], int totalParticles, double boundaries[6])
{
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
        return 1;
    }

    return 0;
}

void cleanup(double (*position)[DIM], double (*velocity)[DIM], double (*acceleration)[DIM],
             double (*oldPosition)[DIM])
{
    delete [] position;
    delete [] velocity;
    delete [] acceleration;
    delete [] oldPosition;
}
