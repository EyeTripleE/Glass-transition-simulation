//double determineVector(const double position1, const double position2, const double higherBound, const double lowerBound)
//{
	//This will need to be modified to accomidate moving borders
	/*if (((position1 - lowerBound >= 4.5f) && (position1 <= higherBound - 4.5f))
		|| (position1 - lowerBound < 4.5f && position2 < (higherBound - (4.5f - position1)))
		|| ((position1 > higherBound - 4.5f) && (position2 >(4.5f - (higherBound - position1)))))
	{
		return (position1 - position2);
	}	
	else if (position1 - lowerBound < 4.5f && position2 > (higherBound - (4.5f - position1)))
	{
		return ((position1 - lowerBound) - (position2 - higherBound));
	}
	else
	{
		return ((position1 - lowerBound) - (position2 + higherBound));
	}*/
	//double ds = position1 - position2;
	//double size = higherBound - lowerBound;
	//return (fabs(ds) > size / 2) ? (ds - copysign(size, ds)) : ds;
//}

double determineVector(const double ds, const double size)
{
	return (fabs(ds) > size / 2) ? (ds - copysign(size, ds)) : ds;
}

double getSigma(const unsigned particlesType1, const unsigned row, const unsigned col)
{
	if (row < particlesType1 && col < particlesType1)
		return 1;
	else if (row >= particlesType1 && col >= particlesType1)
		return 1.4;
	return 1.2;
}
