kernel void calc_acceleration(global double2* position, global double2* acceleration, global double2* wallForce, const double4 boundaries, global double2* energies, global double* virialPotEnergy, global double* wallPE)
{	
	int col = 0;
	int row = get_global_id(0);
	bool val;
	double virialPotEnergyTemp = 0;
	double2 energy = 0;
	double2 forceSums = 0;
	double2 tmpRowVal = position[row];
	double2 size = {boundaries.y - boundaries.x, boundaries.w - boundaries.z};
	double2 vectors;
	double2 virialVectors;
	double2 force;
	double pythagPow6;
	double pythagPow4;
	double pythagPow3;
	double sigmaPow12;
	double sigmaPow6;
	double sigma;
	double pythagorean;
	double instEng;
	/*local double2 positionLocal[2048];
	size_t localID = get_local_id(0);
	size_t numIterations;
	size_t i = 0;
	size_t j;

	while(i < GLOBALSIZE)
	{				
		j = 0;
		//Copy next 2048 elements or all remaining elements into local memory
		while(j < 2048)
		{
			col = i + j + localID;			
			positionLocal[j + localID] = (col < GLOBALSIZE) ? position[col] : 0;
			j += LOCALSIZE;
		}
		barrier(CLK_GLOBAL_MEM_FENCE);

		j = 0;
		numIterations = (i + 2048) > GLOBALSIZE ? (GLOBALSIZE - i) : 2048; 
		while(j < numIterations)
		{
			col = i + j;
			//For if we need to do just one		
			//vectors.x = tmpRowVal.x - positionLocal[j].x;
			//vectors.y = determineVector(tmpRowVal.y - positionLocal[j].y, size.y);

			//For if we need to do both
			vectors = tmpRowVal - positionLocal[j];	
			vectors = (fabs(vectors) > size * 0.5) ? (vectors - copysign(size, vectors)) : vectors;	
					
			pythagorean = dot(vectors, vectors);

			val = (pythagorean < 16 && row != col);
			sigma = getSigma(TYPEONE, row, col);
			sigmaPow6 = pown(sigma, 6);//sigmaPow2*sigmaPow2*sigmaPow2;
			sigmaPow12 = sigmaPow6*sigmaPow6;
			pythagPow3 = pythagorean*pythagorean*pythagorean;
			pythagPow4 = pythagPow3*pythagorean;
			pythagPow6 = pythagPow3*pythagPow3;
			//forceSums += ((24 * vectors * sigmaPow6) / pythagPow4) * ((2 * sigmaPow6 / pythagPow3) - 1);
			force = val ? ((48 * vectors * sigmaPow12) / (pythagPow4 * pythagPow3)) : 0;
			forceSums += force;
			//Multiplying by 2 instead of by four since the potential energy is calculated twice for each pair of particles
			//energy.x += 2 * ((sigmaPow12 / pythagPow6) - (sigmaPow6/pythagPow3));
			energy.x += val ? 2 * (sigmaPow12 / pythagPow6) : 0;	
			virial += val ? dot(vectors, force) : 0;			
			++j;
		}
		i += 2048;
	}*/

	while(col < GLOBALSIZE)
	{
		virialVectors = vectors = tmpRowVal - position[col];//Global memory access
		vectors = (fabs(vectors) > size * 0.5) ? (vectors - copysign(size, vectors)) : vectors;		
		//vectors.y = (fabs(vectors.y) > size.y * 0.5) ? (vectors.y - copysign(size.y, vectors.y)) : vectors.y;		
			
		pythagorean = dot(vectors, vectors);
		val = (pythagorean < 16 && row != col);
		sigma = (row < TYPEONE && col < TYPEONE) ? 1 : ((row >= TYPEONE && col >= TYPEONE) ? 1.4 : 1.2);
		sigmaPow6 = pown(sigma, 6);//sigmaPow2*sigmaPow2*sigmaPow2;
		sigmaPow12 = sigmaPow6*sigmaPow6;
		pythagPow3 = pythagorean*pythagorean*pythagorean;
		pythagPow4 = pythagPow3*pythagorean;
		pythagPow6 = pythagPow3*pythagPow3;
		//forceSums += ((24 * vectors * sigmaPow6) / pythagPow4) * ((2 * sigmaPow6 / pythagPow3) - 1);
		force = val ? ((48 * vectors * sigmaPow12) / (pythagPow4 * pythagPow3)) : 0;
		forceSums += force;
		//Multiplying by 2 instead of by four since the potential energy is calculated twice for each pair of particles
		//energy.x += 2 * ((sigmaPow12 / pythagPow6) - (sigmaPow6/pythagPow3));
		instEng = val ? 2 * (sigmaPow12 / pythagPow6) : 0;
		energy.x += instEng;
		virialPotEnergyTemp += (virialVectors.x == vectors.x) && (virialVectors.y == vectors.y) ? instEng : 0;	
		//virial += val ? dot(vectors, force) : 0;			
		++col;
	}
	
	double privWallPE = 0;
	double2 privWallForce = 0;	

	/*CODE FOR MOVING WALLS*/
	/*double edgeDistancePow11;
	double edgeDistance = tmpRowVal.x - boundaries.x;

	//Within 4 of left bound
	//val = edgeDistance < 4;	
	sigma = row < TYPEONE ? 1 : 1.2;
	sigmaPow12 = powr(sigma, 12);

	//edgeDistancePow11 = val ? powr(edgeDistance, -11) : 0;
	edgeDistancePow11 = powr(edgeDistance, -11);
	//need negative for force on wall 
	privWallForce.x -= 34.01755795*edgeDistancePow11*edgeDistance*sigmaPow12;
	//Subtracting a negative is the same as adding
	forceSums.x -= privWallForce.x; //force goes ===> 
	//virial += edgeDistance*-privWallForce.x;
	privWallPE += 3.092505268*edgeDistancePow11*sigmaPow12;
	
	//Within 4 of right bound... theoretically the sides could get close enough together that they both affect the particle
	edgeDistance = boundaries.y - tmpRowVal.x;
	//val = edgeDistance < 4;
	//edgeDistancePow11 = val ? powr(edgeDistance, -11) : 0;
	edgeDistancePow11 = powr(edgeDistance, -11);
	privWallForce.y += 34.01755795*edgeDistancePow11*edgeDistance*sigmaPow12;
	forceSums.x -= privWallForce.y; //force goes <===
	//virial += edgeDistance*privWallForce.y;
	privWallPE += 3.092505268*edgeDistancePow11*sigmaPow12;
	/*END CODE FOR MOVING WALLS*/

	wallForce[row] = privWallForce;
	acceleration[row] = (row < TYPEONE) ? forceSums : 0.5*forceSums;
	energies[row] = energy;
	virialPotEnergy[row] = virialPotEnergyTemp;
	//innerVirial[row] = virial;
	wallPE[row] = privWallPE;
}

kernel void apply_cooling_boundary(global double2* position, global double2* oldPosition, const double4 boundaries)
{
	size_t index = get_global_id(0);
	double v;
	double vPrime;
	double t;
	double m;
	double2 privatePos = position[index];
	double2 privateOldPos = oldPosition[index];

	//Right
	if(privatePos.x > boundaries.y)
	{
		//Get mass
		m = (index < TYPEONE) ? 1 : 2;
		//Calculate current velocity
		v = (privatePos.x - privateOldPos.x)/TIMESTEP;
		//Calculate time at which particle hits wall
		t = (boundaries.y - privateOldPos.x)/v;
		//Calculate new velocity
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		//Update current position
		privatePos.x = boundaries.y - vPrime*(TIMESTEP - t);
		//Update old position
		privateOldPos.x = boundaries.y + vPrime*t;
	}
	//Left
	else if(privatePos.x < boundaries.x)
	{
		m = (index < TYPEONE) ? 1 : 2;
		v = (privatePos.x - privateOldPos.x)/TIMESTEP;
		t = (boundaries.x - privateOldPos.x)/v;
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		privatePos.x = boundaries.x + vPrime*(TIMESTEP - t);
		privateOldPos.x = boundaries.x - vPrime*t;
	}
	//Top
	if(privatePos.y > boundaries.w)
	{
		m = (index < TYPEONE) ? 1 : 2;
		v = (privatePos.y - privateOldPos.y)/TIMESTEP;
		t = (boundaries.w - privateOldPos.y)/v;
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		privatePos.y = boundaries.w - vPrime*(TIMESTEP - t);
		privateOldPos.y = boundaries.w + vPrime*t;
	}
	//Bottom
	else if(privatePos.y < 0)
	{
		m = (index < TYPEONE) ? 1 : 2;
		v = (privatePos.y - privateOldPos.y)/TIMESTEP;
		t = (-privateOldPos.y)/v;
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		privatePos.y = vPrime*(TIMESTEP - t);
		privateOldPos.y = -vPrime*t;
	}

	position[index] = privatePos;
	oldPosition[index] = privateOldPos;
}

kernel void apply_cooling_membrane(global double2* position, global double2* oldPosition, const double4 boundaries)
{
	size_t index = get_global_id(0);
	double v;
	double vPrime;
	double t;
	double m;
	double2 privatePos = position[index];
	double2 privateOldPos = oldPosition[index];

	if(privatePos.x > boundaries.y)
	{
		m = (index < TYPEONE) ? 1 : 2;
		v = (privatePos.x - privateOldPos.x)/TIMESTEP;
		t = (boundaries.y - privateOldPos.x)/v;
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		privatePos.x = boundaries.x + vPrime*(TIMESTEP - t);
		privateOldPos.x = boundaries.x - vPrime*t;
	}
	else if(privatePos.x < boundaries.x)
	{
		m = (index < TYPEONE) ? 1 : 2;
		v = (privatePos.x - privateOldPos.x)/TIMESTEP;
		t = (boundaries.x - privateOldPos.x)/v;
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		privatePos.x = boundaries.y - vPrime*(TIMESTEP - t);
		privateOldPos.x = boundaries.y + vPrime*t;
	}

	if(privatePos.y > boundaries.w)
	{
		m = (index < TYPEONE) ? 1 : 2;
		v = (privatePos.y - privateOldPos.y)/TIMESTEP;
		t = (boundaries.w - privateOldPos.y)/v;
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		privatePos.y = vPrime*(TIMESTEP - t);
		privateOldPos.y = -vPrime*t;
	}
	else if(privatePos.y < 0)
	{
		m = (index < TYPEONE) ? 1 : 2;
		v = (privatePos.y - privateOldPos.y)/TIMESTEP;
		t = (-privateOldPos.y)/v;
		vPrime = sqrt((-1/(BETADIV2*m))*log1p(-exp(-v*v*BETADIV2*m)));
		privatePos.y = boundaries.w - vPrime*(TIMESTEP - t);
		privateOldPos.y = boundaries.w + vPrime*t;
	}

	position[index] = privatePos;
	oldPosition[index] = privateOldPos;
}

kernel void apply_solid_boundary(global double2* position, global double2* oldPosition, const double4 boundaries)
{
	size_t index = get_global_id(0);
	double ds;
	double2 privatePos = position[index];
	double2 privateOldPos = oldPosition[index];

	if(privatePos.x > boundaries.y)
	{
		ds = privatePos.x - boundaries.y;
		privatePos.x -= (ds + ds);
		ds = boundaries.y - privateOldPos.x;
		privateOldPos.x += (ds + ds);
	}
	else if(privatePos.x < boundaries.x)
	{
		ds = boundaries.x - privatePos.x;
		privatePos.x += (ds + ds);
		ds = privateOldPos.x - boundaries.x;
		privateOldPos.x -= (ds + ds);
	}

	if(privatePos.y > boundaries.w)
	{
		ds = privatePos.y - boundaries.w;
		privatePos.y -= (ds + ds);
		ds = boundaries.w - privateOldPos.y;
		privateOldPos.y += (ds + ds);
	}
	else if(privatePos.y < 0)
	{
		ds = -privatePos.y;
		privatePos.y += (ds + ds);
		ds = privateOldPos.y;
		privateOldPos.y -= (ds + ds);
	}

	position[index] = privatePos;
	oldPosition[index] = privateOldPos;
}

kernel void apply_periodic_boundary(global double2* position, global double2* oldPosition, const double4 boundaries)
{
	size_t index = get_global_id(0);
	double holder;
	double intervals = boundaries.y - boundaries.x;

	//Copy from global memory
	double2 privatePos = position[index];
	double2 privateOldPos = oldPosition[index];

	//Need this in case particle jumps over piston
	//if (privatePos.x > boundaries.y || privatePos.x < boundaries.x)
	//{
		//holder = (privatePos.x > boundaries.y || privatePos.x < boundaries.x) ? copysign(intervals, privatePos.x - intervals) : 0;
		//privatePos.x -= holder;
		//privateOldPos.x -= holder;
	//}

	/*while (privatePos.x > boundaries.y)
	{
		privatePos.x -= (intervals);
		privateOldPos.x -= (intervals);
	}
	while (privatePos.x < boundaries.x)
	{
		privatePos.x += (intervals);
		privateOldPos.x += (intervals);
	}*/

	//if (privatePos.y > boundaries.w || privatePos.y < 0)
	//{
		holder = (privatePos.x > boundaries.y || privatePos.x < boundaries.x) ? copysign(intervals, privatePos.x - intervals) : 0;
		privatePos.x -= holder;
		privateOldPos.x -= holder;

		holder = (privatePos.y > boundaries.w || privatePos.y < 0) ? copysign(boundaries.w, privatePos.y) : 0;
		privatePos.y -= holder;
		privateOldPos.y -= holder;
	//}

	/*while (privatePos.y > boundaries.w)
	{
		privatePos.y -= boundaries.w;
		privateOldPos.y -= boundaries.w;
	}
	while (privatePos.y < 0)
	{
		privatePos.y += boundaries.w;
		privateOldPos.y += boundaries.w;
	}*/
	position[index] = privatePos;
	oldPosition[index] = privateOldPos;
}

kernel void euler_kernel(global double2* position, global double2* velocity, global double2* acceleration, global double2* oldPosition)
{
	size_t index = get_global_id(0);
	//Copy to private memory
	double2 tmpPos = position[index];
	double2 tmpAcc = acceleration[index];
	double2 tmpVel = velocity[index];
	//Execute Euler
	oldPosition[index] = tmpPos;
	//velocity[index] = tmpVel + timestep*tmpAcc;
	position[index] = tmpPos + tmpVel * TIMESTEP + 0.5*TIMESTEP*TIMESTEP*tmpAcc;
}

kernel void verlet_kernel(global double2* position, global double2* acceleration, global double2* oldPosition, global double2* energies)
{
	//Get the index of the work-item
	size_t index = get_global_id(0);

	//Vector Verlet Method
	double2 curPos = position[index];
	double2 oldPos = oldPosition[index];
	double2 energy = 0;//{ 0, 0 };

	//find displacement
	double2 displacement = curPos - oldPos;
	//find future displacement
	displacement += TIMESTEP * TIMESTEP * acceleration[index];
	//update oldPosition
	double2 oldoldPos = oldPos;
	oldPos = curPos;  
	//find future position
	curPos += displacement;
	position[index] = curPos;
	//determine velocity at the previous timestep
	//Subtractive cancellation here potentially
	double2 vel = ((curPos - oldoldPos) / (2 * TIMESTEP));
	oldPosition[index] = oldPos;

	energy.y = (index < TYPEONE) ? 0.5*dot(vel, vel) : dot(vel, vel);
	//Need += otherwise will wipe out the potential energy
	energies[index] += energy;
}

kernel void reduceDouble2(const global double2* inputData, global double2* outputData, local double2* tempData)
{
	size_t global_id = get_global_id(0);
	int local_id = get_local_id(0);
	int i = 0;
	double2 accumulator = 0;//{0, 0};

	while(i < STEPSIZE)
	{
		//Global memory access... inefficient
		//accumulator += inputData[local_id*STEPSIZE + i];
		accumulator += inputData[mad24(local_id, STEPSIZE, i)];
		++i;
	}
	tempData[local_id] = accumulator;	

	//We should be left only with an array of length local_size which we can parallel reduce	

	//The local size must be a power of two and the size of inputData must be a multiple of this.
	i = LOCALSIZE >> 1;
	while(i > 0)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		tempData[local_id] += (local_id < i) ? tempData[local_id + i] : 0;
		i >>= 1;
	}
	//Should wind up with tempData[0] being the sum of the local chunk
	
	if (global_id == 0)
	{
		//Copy to the output
		outputData[0] = tempData[0];
	}
}

kernel void reduceDouble(const global double* inputData, global double* outputData, local double* tempData)
{
	size_t global_id = get_global_id(0);
	int local_id = get_local_id(0);
	int i = 0;
	double accumulator = 0;

	while(i < STEPSIZE)
	{
		//accumulator += inputData[local_id*STEPSIZE + i];
		accumulator += inputData[mad24(local_id, STEPSIZE, i)];
		++i;
	}
	tempData[local_id] = accumulator;

	//We should be left only with an array of length local_size which we can parallel reduce	

	//The local size must be a power of two and the size of inputData must be a multiple of this.
	i = LOCALSIZE >> 1;
	while(i > 0)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		tempData[local_id] += (local_id < i) ? tempData[local_id + i] : 0;
		i >>= 1;
	}
	//Should wind up with tempData[0] being the sum of the local chunk

	//Copy to the output
	if (global_id == 0)
	{
		outputData[0] = tempData[0];
	}
}

kernel void reduceUint(const global uint* inputData, global uint* outputData, local uint* tempData)
{
	size_t global_id = get_global_id(0);
	int local_id = get_local_id(0);
	int i = 0;
	uint accumulator = 0;

	while(i < STEPSIZE)
	{
		//accumulator += inputData[local_id*STEPSIZE + i];
		accumulator += inputData[mad24(local_id, STEPSIZE, i)];
		++i;
	}
	tempData[local_id] = accumulator;

	//We should be left only with an array of length local_size which we can parallel reduce	

	//The local size must be a power of two and the size of inputData must be a multiple of this.
	i = LOCALSIZE >> 1;
	while(i > 0)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		tempData[local_id] += (local_id < i) ? tempData[local_id + i] : 0;
		i >>= 1;
	}
	//Should wind up with tempData[0] being the sum of the local chunk

	//Copy to the output
	if (global_id == 0)
	{
		outputData[0] = tempData[0];
	}
}

kernel void regionDensityKernel(const global double2* position, global uint* result, const double4 check)
{
	size_t index = get_global_id(0);
	double2 privPos = position[index];
	result[index] = ((privPos.x >= check.x) &&
					(privPos.x <= check.y) &&
					(privPos.y >= check.z) &&
					(privPos.y <= check.w)) ? 1 : 0;
}

//Not likely possible to do because of race conditions
/*kernel void pairCorrelationKernel(global uint* bins1to1, global uint* bins1to2, global uint* bins2to2, global double3* adjustedBins, global double2* position, const uint totalParticles, const uint numParticlesType1, const double4 boundaries, const uint maxRange, const double binSize)
{
	//READ_WRITE or WRITE_ONLY?
	//double innerRadius;
	//double outerRadius;
	//double area;
	uint global_id = get_global_id(0);
	//int localID = get_local_id(0);
	//col pretty much functions as 'i'
	uint col;
	double2 vectors;
	double2 tmpRow = position[global_id];
	double vect_length;
	uint binNumber;
	//double density = ((double)totalParticles) / (boundaries.x * boundaries.y);

	uint writeLocation;
	size_t globalsize = get_global_size(0);

	//Find way to zero buffers
	//Buffers can be any length
	for (col = 0; col < maxRange; col += globalsize)
	{
		writeLocation = col + global_id;
		if (writeLocation < maxRange)
		{
			bins1to1[writeLocation] = 0;
			bins1to2[writeLocation] = 0;
			bins2to2[writeLocation] = 0;
		}
	}
	//Make sure global memory all zeroed
	barrier(CLK_GLOBAL_MEM_FENCE);

	for (col = 1; col < totalParticles; ++col)
	{
		//Makes for a triangular 
		if (col < global_id)
		{
			vectors.x = determineVector(tmpRow.x, position[col].x, boundaries.x);
			vectors.y = determineVector(tmpRow.y, position[col].y, boundaries.y);
			vect_length = length(vectors);
			binNumber = (int)(vect_length / binSize);

			if (global_id < numParticlesType1 && col < numParticlesType1)
			{
				atomic_inc(&bins1to1[binNumber]);
			}
			else if (global_id >= numParticlesType1 && col >= numParticlesType1)
			{
				atomic_inc(&bins1to2[binNumber]);
			}
			else
			{
				atomic_inc(&bins2to2[binNumber]);
			}
		}
	}

	barrier(CLK_GLOBAL_MEM_FENCE);

	for (col = 0; col < maxRange; col += globalsize)
	{
		writeLocation = col + global_id;
		if (writeLocation < maxRange)
		{
			adjustedBins[writeLocation].x = bins1to1[writeLocation];
			adjustedBins[writeLocation].y = bins1to2[writeLocation];
			adjustedBins[writeLocation].z = bins2to2[writeLocation];
		}
	}

	//barrier(CLK_GLOBAL_MEM_FENCE);
	
	for (col = 0; col < maxRange; col += globalsize)
	{
		writeLocation = col + global_id;//row is the same as the global ID
		if (writeLocation < maxRange)
		{
			outerRadius = (writeLocation + 1) * binSize;
			outerRadius *= outerRadius;
			innerRadius = (writeLocation) * binSize;
			innerRadius *= innerRadius;
			area = M_PI_F*(outerRadius - innerRadius);
			//area = pi*(((writeLocation + 1) * bsize) * ((writeLocation + 1) * bsize) - (writeLocation * bsize) * (writeLocation * bsize));
			adjustedBins[writeLocation] /= area; //Note: write_location is wrong when assigned to array[writeLocation]
		}
	}
}*/