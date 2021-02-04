#include "SS.h"

double spatial_kernel(double*q, double*p, statistics& stat)
{
	double value;

	if (stat.k_type_s == 0) //Triangular kernel
	{
		value = 1 - stat.gamma_s*fabs(sqrt((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1])));

		if (value < 0)
			return 0;
		else
			return value;
	}

	if (stat.k_type_s == 1) //Epanechnikov kernel
	{
		value = 1 - stat.gamma_s*stat.gamma_s*((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1]));

		if (value < 0)
			return 0;
		else
			return value;
	}

	if (stat.k_type_s == 2) //Quartic kernel
	{
		value = 1 - stat.gamma_s*stat.gamma_s*((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1]));

		if (value < 0)
			return 0;
		else
			return value * value;
	}

	return -inf;
}

double temporal_kernel(double*q, double*p, statistics& stat)
{
	double value;
	
	if (stat.k_type_t == 0) //Triangular kernel
	{
		value = 1 - stat.gamma_t*fabs(q[2] - p[2]);
		if (value < 0)
			return 0;
		else
			return value;
	}

	if (stat.k_type_t == 1) //Epanechnikov kernel
	{
		value = 1 - stat.gamma_t*stat.gamma_t*(q[2] - p[2])*(q[2] - p[2]);
		if (value < 0)
			return 0;
		else
			return value;
	}

	if (stat.k_type_t == 2) //Quartic kernel
	{
		value = 1 - stat.gamma_t*stat.gamma_t*(q[2] - p[2])*(q[2] - p[2]);
		if (value < 0)
			return 0;
		else
			return value * value;
	}

	return -inf;
}

double SCAN(statistics& stat)
{
	double incr_value;

	incr_value = 0;
	for (int i = 0; i < stat.n; i++)
		incr_value += spatial_kernel(stat.q, stat.featureVector[i], stat)*temporal_kernel(stat.q, stat.featureVector[i], stat);

	return incr_value;
}

double SCAN_RQS_result(statistics& stat, Tree& tree)
{
	double incr_value;
	int id;

	incr_value = 0;
	for (int i = 0; i < (int)tree.range_result_idList.size(); i++)
	{
		id = tree.range_result_idList[i];
		incr_value += spatial_kernel(stat.q, stat.featureVector[id], stat)*temporal_kernel(stat.q, stat.featureVector[id], stat);
	}

	return incr_value;
}