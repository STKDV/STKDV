#pragma once
#ifndef INIT_VISUAL_H
#define INIT_VISUAL_H

#include "Library.h"

const double inf = 999999999999;
const double small_epsilon = 0.000000000000001;
const double pi = 3.14159265358979323846;

struct statistics
{
	int n; //number of points in datasets
	int n_x; //number of discrete region in the x-axis, e.g., 100
	int n_y; //number of discrete region in the y-axis, e.g., 100
	int n_t; //number of discrete region in the t-axis, e.g., 100
	double x_L, x_U; //lower and upper region of data points in x-axis
	double y_L, y_U; //lower and upper region of data points in y-axis
	double t_L, t_U; //lower and upper region of data points in t-axis
	double incr_x; //incremental value in the x-axis
	double incr_y; //incremental value in the y-axis
	double incr_t; //incremental value in the t-axis
	const int dim = 3; //dim = 3 (3d visualization)
	int method; //chosen method
	int k_type_s; //kernel type 0: Gaussian kernel, 1: Epanechnikov kernel
	int k_type_t; //kernel type 0: Gaussian kernel, 1: Epanechnikov kernel
	double gamma_ratio_s; //ratio = 1 means using Scott's rule for spatial kernel
	double gamma_ratio_t; //ratio = 1 means using Scott's rule for temporal kernel
	double gamma_s; //Control the bandwidth of spatial kernel
	double gamma_t; //Control the bandwidth of temporal kernel
	char*dataFileName;
	char*outputFileName;
	double b; //From Gan17 (SIGMOD17) 

	double*q;
	double**featureVector; //feature vector of all data points (d = 0, 1 and 2 represent longitude, latitude and timestamp respectively)
	double***out_tensor; //output visualization

	//Used in RQS method
	double**query_boundary;

	//Used in SWS method
	double*sliding_window;
	double**sorted_featureVector;
	double*sliding_window_L; //Used in triangular kernel
	double*sliding_window_R; //Used in triangular kernel

	//Used in progressive visualization framework
	bool is_progressive;
	char*prev_level_resultfileName;
};

//distance computation
double s_dist(double*q, double*p);
double t_dist(double*q, double*p);
double e_dist(double*q, double*p);

//offline
void initTensor(statistics& stat);
void initStat(int argc, char**argv, statistics& stat);
void extract_FeatureVector(statistics& stat);
void updateRegion(statistics& stat);
void update_incr_value(statistics& stat);
void init_visual(int argc, char**argv, statistics& stat);
void obtain_bandwidth_values(statistics& stat);
void load_result_file(statistics& stat); //Used in progressive visualization framework

//online
void obtain_q(int x_index, int y_index, int t_index, statistics& stat);

//output result
void save_result_to_file(statistics& stat);

#endif