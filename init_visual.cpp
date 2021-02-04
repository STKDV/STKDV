#include "init_visual.h"

double s_dist(double*q, double*p)
{
	double euclid_dist;
	euclid_dist = sqrt((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1]));

	return euclid_dist;
}

double t_dist(double*q, double*p)
{
	double euclid_dist;
	euclid_dist = fabs(q[2] - p[2]);
	return euclid_dist;
}

double e_dist(double*q, double*p)
{
	double euclid_dist;
	euclid_dist = sqrt((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1]) + (q[2] - p[2])*(q[2] - p[2]));
	return euclid_dist;
}

void initTensor(statistics& stat)
{
	stat.out_tensor = new double**[stat.n_x];
	for (int x = 0; x < stat.n_x; x++)
	{
		stat.out_tensor[x] = new double*[stat.n_y];
		for (int y = 0; y < stat.n_y; y++)
			stat.out_tensor[x][y] = new double[stat.n_t];
	}
}

void initStat(int argc, char**argv, statistics& stat)
{
	stat.dataFileName = argv[1];
	stat.outputFileName = argv[2];
	stat.method = atoi(argv[3]);
	stat.n_x = atoi(argv[4]);
	stat.n_y = atoi(argv[5]);
	stat.n_t = atoi(argv[6]);
	stat.k_type_s = atoi(argv[7]);
	stat.k_type_t = atoi(argv[8]);
	stat.gamma_ratio_s = atof(argv[9]);
	stat.gamma_ratio_t = atof(argv[10]);
	stat.b = atof(argv[11]);
	stat.is_progressive = atoi(argv[12]);

	if (stat.is_progressive == true)
		stat.prev_level_resultfileName = argv[13];

	//debug
	/*stat.dataFileName = (char*)"../../../Datasets/Testing/Testing";
	stat.outputFileName = (char*)"./Results/Testing_M4";
	stat.method = 4;
	stat.n_x = 16;
	stat.n_y = 16;
	stat.n_t = 16;
	stat.k_type_s = 1;
	stat.k_type_t = 1;
	stat.gamma_ratio_s = 1;
	stat.gamma_ratio_t = 1;*/
}

void extract_FeatureVector(statistics& stat)
{
	//load data to feature array
	fstream file;
	string temp_string;

	file.open(stat.dataFileName);
	if (file.is_open() == false)
	{
		cout << "Cannot Open File!" << endl;
		exit(1);
	}

	file >> temp_string;
	file >> stat.n;
	
	stat.featureVector = new double*[stat.n];
	for (int i = 0; i < stat.n; i++)
		stat.featureVector[i] = new double[stat.dim];

	file >> temp_string; file >> temp_string; file >> temp_string;
	for (int i = 0; i < stat.n; i++)
		for (int d = 0; d < stat.dim; d++)
			file >> stat.featureVector[i][d];

	file.close();
}

void updateRegion(statistics& stat)
{
	stat.x_L = inf; stat.x_U = -inf;
	stat.y_L = inf;	stat.y_U = -inf;
	stat.t_L = inf;	stat.t_U = -inf;

	for (int i = 0; i < stat.n; i++)
	{
		if (stat.featureVector[i][0] < stat.x_L)
			stat.x_L = stat.featureVector[i][0];
		if (stat.featureVector[i][0] > stat.x_U)
			stat.x_U = stat.featureVector[i][0];

		if (stat.featureVector[i][1] < stat.y_L)
			stat.y_L = stat.featureVector[i][1];
		if (stat.featureVector[i][1] > stat.y_U)
			stat.y_U = stat.featureVector[i][1];

		if (stat.featureVector[i][2] < stat.t_L)
			stat.t_L = stat.featureVector[i][2];
		if (stat.featureVector[i][2] > stat.t_U)
			stat.t_U = stat.featureVector[i][2];
	}
}

void update_incr_value(statistics& stat)
{
	stat.q = new double[stat.dim];
	stat.incr_x = (stat.x_U - stat.x_L) / stat.n_x;
	stat.incr_y = (stat.y_U - stat.y_L) / stat.n_y;
	stat.incr_t = (stat.t_U - stat.t_L) / stat.n_t;
}

void init_visual(int argc, char**argv, statistics& stat)
{
	initStat(argc, argv, stat);
	initTensor(stat);
	extract_FeatureVector(stat);
	updateRegion(stat);
	update_incr_value(stat);
	obtain_bandwidth_values(stat);

	if (stat.is_progressive == true)
		load_result_file(stat);
}

void save_result_to_file(statistics& stat)
{
	fstream outputFile;
	outputFile.open(stat.outputFileName, ios::in | ios::out | ios::trunc);

	outputFile << "x_L " << stat.x_L << endl;
	outputFile << "x_U " << stat.x_U << endl;
	outputFile << "y_L " << stat.y_L << endl;
	outputFile << "y_U " << stat.y_U << endl;
	outputFile << "t_L " << stat.t_L << endl;
	outputFile << "t_U " << stat.t_U << endl;
	outputFile << "n_x " << stat.n_x << endl;
	outputFile << "n_y " << stat.n_y << endl;
	outputFile << "n_t " << stat.n_t << endl;

	for (int x = 0; x < stat.n_x; x++)
		for (int y = 0; y < stat.n_y; y++)
			for (int t = 0; t < stat.n_t; t++)
				outputFile << stat.out_tensor[x][y][t] << endl;

	outputFile.close();
}

void load_result_file(statistics& stat)
{
	fstream outputFile;
	int temp;

	outputFile.open(stat.outputFileName);

	outputFile >> temp; outputFile >> temp; outputFile >> temp;
	outputFile >> temp; outputFile >> temp; outputFile >> temp;
	outputFile >> temp; outputFile >> temp; outputFile >> temp;

	for (int x = 0; x < stat.n_x; x++)
		for (int y = 0; y < stat.n_y; y++)
			for (int t = 0; t < stat.n_t; t++)
				outputFile >> stat.out_tensor[x][y][t];

	outputFile.close();
}

void obtain_q(int x_index, int y_index, int t_index, statistics& stat)
{
	stat.q[0] = stat.x_L + x_index * stat.incr_x;
	stat.q[1] = stat.y_L + y_index * stat.incr_y;
	stat.q[2] = stat.t_L + t_index * stat.incr_t;
}

void obtain_bandwidth_values(statistics& stat)
{
	double sum_x, sum_y, sum_t;
	double mean_x, mean_y, mean_t;
	double sd_x, sd_y, sd_t;
	double h_x, h_y, h_t;

	//Using Scott's rule to obtain the bandwidth for both spatial and temporal kernels

	//Compute mean value
	sum_x = 0; sum_y = 0; sum_t = 0;
	for (int i = 0; i < stat.n; i++)
	{
		sum_x += stat.featureVector[i][0];
		sum_y += stat.featureVector[i][1];
		sum_t += stat.featureVector[i][2];
	}
	mean_x = sum_x / stat.n; mean_y = sum_y / stat.n; mean_t = sum_t / stat.n;

	sd_x = 0; sd_y = 0; sd_t = 0;
	for (int i = 0; i < stat.n; i++)
	{
		sd_x += (stat.featureVector[i][0] - mean_x)*(stat.featureVector[i][0] - mean_x) / (stat.n - 1);
		sd_y += (stat.featureVector[i][1] - mean_y)*(stat.featureVector[i][1] - mean_y) / (stat.n - 1);
		sd_t += (stat.featureVector[i][2] - mean_t)*(stat.featureVector[i][2] - mean_t) / (stat.n - 1);
	}
	sd_x = sqrt(sd_x); sd_y = sqrt(sd_y); sd_t = sqrt(sd_t);
	h_x = stat.b * pow((double)stat.n, -1.0 / 6.0)*sd_x; //d = 2
	h_y = stat.b * pow((double)stat.n, -1.0 / 6.0)*sd_y; //d = 2
	h_t = stat.b * pow((double)stat.n, -1.0 / 5.0)*sd_t; //d = 1

	stat.gamma_s = stat.gamma_ratio_s / (sqrt(h_x*h_x + h_y * h_y));
	stat.gamma_t = stat.gamma_ratio_t / h_t;
}