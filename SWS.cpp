#include "SWS.h"

void init_SWS(statistics& stat)
{
	vector<index_time_pair> pair_vector;
	index_time_pair pair;
	int cur_index;

	stat.sorted_featureVector = new double*[stat.n];
	for (int i = 0; i < stat.n; i++)
	{
		pair.index = i;
		pair.time = stat.featureVector[i][2];
		pair_vector.push_back(pair);
		stat.sorted_featureVector[i] = new double[stat.dim];
	}

	sort(pair_vector.begin(), pair_vector.end());

	for (int i = 0; i < stat.n; i++)
	{
		cur_index = pair_vector[i].index;
		for (int d = 0; d < stat.dim; d++)
			stat.sorted_featureVector[i][d] = stat.featureVector[cur_index][d];
	}

	if (stat.k_type_t == 0) //Triangular kernel
	{
		stat.sliding_window_L = new double[2];
		stat.sliding_window_R = new double[2];
	}
	if (stat.k_type_t == 1) //Epanechnikov kernel
		stat.sliding_window = new double[3];
	if (stat.k_type_t == 2) //Quartic kernel
		stat.sliding_window = new double[5];
}

double compute_init_window_density(statistics& stat, win_status& win)
{
	win.start_t_win_val = stat.q[2] - win.t_bandwidth;
	win.end_t_win_val = stat.q[2] + win.t_bandwidth;
	win.start_t_win_pos = 0;
	double win_density_value;
	double kernel_s_value;
	double t_p_pow;
	bool isStart = false;
	bool isEnd = false;
	//Used in Quartic kernel
	double gamma_t_pow_2, gamma_t_pow_4;
	double t_q_pow_2, t_q_pow_3, t_q_pow_4;

	//initialization
	if (stat.k_type_t == 1) //Epanechnikov kernel
	{
		for (int w = 0; w < 3; w++)
			stat.sliding_window[w] = 0;
	}
	if (stat.k_type_t == 2) //Quartic kernel
	{
		for (int w = 0; w < 5; w++)
			stat.sliding_window[w] = 0;
	}

	for (int i = 0; i < stat.n; i++)
	{
		if (isStart == false)
		{
			if (stat.sorted_featureVector[i][2] > win.start_t_win_val) //This is the first point that is bigger than the starting point of the interval
			{
				isStart = true;
				win.start_t_win_pos = i;
			}
		}

		if (isEnd == false)
		{
			if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
			{
				isEnd = true;
				win.end_t_win_pos = i - 1;
			}
			else
			{
				if (isStart == true) //This point is inside the interval
				{
					kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[i], stat);
					t_p_pow = 1;

					stat.sliding_window[0] += kernel_s_value;
					if (stat.k_type_t == 1) //Epanechnikov kernel
					{
						for (int w = 1; w < 3; w++)
						{
							t_p_pow *= stat.sorted_featureVector[i][2];
							stat.sliding_window[w] += t_p_pow * kernel_s_value;
						}
					}
					if (stat.k_type_t == 2) //Quartic kernel
					{
						for (int w = 1; w < 5; w++)
						{
							t_p_pow *= stat.sorted_featureVector[i][2];
							stat.sliding_window[w] += t_p_pow * kernel_s_value;
						}
					}
				}

				if (i == stat.n - 1) //The last point
					win.end_t_win_pos = stat.n - 1;
			}
		}

		if (isStart == true && isEnd == true)
			break;
	}

	if (stat.k_type_t == 1) //Epanechnikov kernel
		win_density_value = (1 - stat.gamma_t*stat.gamma_t*stat.q[2] * stat.q[2])*stat.sliding_window[0]
		+ 2 * stat.gamma_t*stat.gamma_t*stat.q[2] * stat.sliding_window[1]
		- stat.gamma_t*stat.gamma_t*stat.sliding_window[2];

	if (stat.k_type_t == 2) //Quartic kernel
	{
		gamma_t_pow_2 = stat.gamma_t*stat.gamma_t; 
		gamma_t_pow_4 = gamma_t_pow_2 * gamma_t_pow_2;
		t_q_pow_2 = stat.q[2] * stat.q[2];
		t_q_pow_3 = t_q_pow_2 * stat.q[2];
		t_q_pow_4 = t_q_pow_3 * stat.q[2];
		win_density_value = (1 - 2 * gamma_t_pow_2*t_q_pow_2 + gamma_t_pow_4 * t_q_pow_4)*stat.sliding_window[0]
			+ (4 * gamma_t_pow_2*stat.q[2] - 4 * gamma_t_pow_4*t_q_pow_3)*stat.sliding_window[1]
			+ (6 * gamma_t_pow_4*t_q_pow_2 - 2 * gamma_t_pow_2)*stat.sliding_window[2]
			- 4 * gamma_t_pow_4*stat.q[2] * stat.sliding_window[3] + gamma_t_pow_4 * stat.sliding_window[4];
	}

	return win_density_value;
}

void update_sliding_window(statistics& stat, vector<int>& index_set, bool is_positive)
{
	double weight;
	double kernel_s_value;
	double t_p_pow;
	int id;

	if (is_positive == true)
		weight = 1;
	else
		weight = -1;

	for (int i = 0; i < (int)index_set.size(); i++)
	{
		id = index_set[i];

		kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);
		stat.sliding_window[0] += weight * kernel_s_value;

		if (stat.k_type_t == 1) //Epanechnikov kernel
		{
			t_p_pow = 1;
			for (int w = 1; w < 3; w++)
			{
				t_p_pow *= stat.sorted_featureVector[id][2];
				stat.sliding_window[w] += weight * t_p_pow * kernel_s_value;
			}
		}

		if (stat.k_type_t == 2) //Quartic kernel
		{
			t_p_pow = 1;
			for (int w = 1; w < 5; w++)
			{
				t_p_pow *= stat.sorted_featureVector[id][2];
				stat.sliding_window[w] += weight * t_p_pow * kernel_s_value;
			}
		}
	}
}

double incr_update_window_density(statistics& stat, win_status& win)
{
	double win_density_value;
	vector<int> D_set;
	vector<int> I_set;
	int cur_index;
	bool isStart = false;
	bool isEnd = false;
	//Used in Quartic kernel
	double gamma_t_pow_2, gamma_t_pow_4;
	double t_q_pow_2, t_q_pow_3, t_q_pow_4;

	win.start_t_win_val_prev = win.start_t_win_val;
	win.end_t_win_val_prev = win.end_t_win_val;
	win.start_t_win_val = stat.q[2] - win.t_bandwidth;
	win.end_t_win_val = stat.q[2] + win.t_bandwidth;

	cur_index = win.start_t_win_pos;
	for (int i = cur_index; i < stat.n; i++)
	{
		if (isStart == false)
		{
			if (stat.sorted_featureVector[i][2] > win.start_t_win_val)
			{
				win.start_t_win_pos = i;
				isStart = true;
			}
		}

		if (isStart == true)
			break;

		if (stat.sorted_featureVector[i][2] <= min(win.end_t_win_val_prev, win.start_t_win_val))
			D_set.push_back(i);
	}

	cur_index = win.end_t_win_pos;
	for (int i = cur_index; i < stat.n; i++)
	{
		if (isEnd == false)
		{
			if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
			{
				win.end_t_win_pos = i - 1;
				isEnd = true;
			}
		}

		if (isEnd == true)
			break;

		if (stat.sorted_featureVector[i][2] > max(win.end_t_win_val_prev, win.start_t_win_val))
			I_set.push_back(i);
	}

	//update sliding window (Equation 9)
	update_sliding_window(stat, D_set, false);
	update_sliding_window(stat, I_set, true);

	if (stat.k_type_t == 1) //Epanechnikov kernel
		win_density_value = (1 - stat.gamma_t*stat.gamma_t*stat.q[2] * stat.q[2])*stat.sliding_window[0]
		+ 2 * stat.gamma_t*stat.gamma_t*stat.q[2] * stat.sliding_window[1]
		- stat.gamma_t*stat.gamma_t*stat.sliding_window[2];
	if (stat.k_type_t == 2) //Quartic kernel
	{
		gamma_t_pow_2 = stat.gamma_t*stat.gamma_t;
		gamma_t_pow_4 = gamma_t_pow_2 * gamma_t_pow_2;
		t_q_pow_2 = stat.q[2] * stat.q[2];
		t_q_pow_3 = t_q_pow_2 * stat.q[2];
		t_q_pow_4 = t_q_pow_3 * stat.q[2];
		win_density_value = (1 - 2 * gamma_t_pow_2*t_q_pow_2 + gamma_t_pow_4 * t_q_pow_4)*stat.sliding_window[0]
			+ (4 * gamma_t_pow_2*stat.q[2] - 4 * gamma_t_pow_4*t_q_pow_3)*stat.sliding_window[1]
			+ (6 * gamma_t_pow_4*t_q_pow_2 - 2 * gamma_t_pow_2)*stat.sliding_window[2]
			- 4 * gamma_t_pow_4*stat.q[2] * stat.sliding_window[3] + gamma_t_pow_4 * stat.sliding_window[4];
	}

	return win_density_value;
}

double compute_init_window_density_tri(statistics& stat, win_status& win)
{
	win.start_t_win_val = stat.q[2] - win.t_bandwidth;
	win.end_t_win_val = stat.q[2] + win.t_bandwidth;
	win.start_t_win_pos = 0;
	double win_density_value;
	double kernel_s_value;
	double t_p_pow;
	bool isStart = false;
	bool isCenter = false;

	for (int w = 0; w < 2; w++)
	{
		stat.sliding_window_L[w] = 0;
		stat.sliding_window_R[w] = 0;
	}

	for (int i = 0; i < stat.n; i++)
	{
		if (isStart == false)
		{
			if (stat.sorted_featureVector[i][2] > win.start_t_win_val) //This is the first point that is bigger than the starting point of the interval
			{
				isStart = true;
				win.start_t_win_pos = i;
			}
		}

		if (isCenter == false && isStart == true)
		{
			if (stat.sorted_featureVector[i][2] > stat.q[2])
			{
				win.center_t_win_pos = i;
				isCenter = true;
			}

			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[i], stat);
			stat.sliding_window_L[0] += kernel_s_value;
			stat.sliding_window_L[1] += stat.sorted_featureVector[i][2] * kernel_s_value;
		}

		if (isCenter == true)
		{
			if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
			{
				win.end_t_win_pos = i - 1;
				break;
			}

			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[i], stat);
			stat.sliding_window_R[0] += kernel_s_value;
			stat.sliding_window_R[1] += stat.sorted_featureVector[i][2] * kernel_s_value;
		}
	}

	win_density_value = (stat.sliding_window_L[0] + stat.sliding_window_R[0])
		- stat.gamma_t*(stat.q[2] * stat.sliding_window_L[0] - stat.sliding_window_L[1]
			+ stat.sliding_window_R[1] - stat.q[2] * stat.sliding_window_R[0]);

	return win_density_value;
}

double incr_update_window_density_tri(statistics& stat, win_status& win)
{
	double win_density_value;
	vector<int> D_set;
	vector<int> I_set;
	vector<int> C_set; //C(t_q, t_{q}_n)
	vector<int> S_set;
	vector<int> A_set;
	int cur_index;
	bool isStart = false;
	bool isCenter = false;
	bool isEnd = false;
	double kernel_s_value;
	int id;

	win.start_t_win_val_prev = win.start_t_win_val;
	win.end_t_win_val_prev = win.end_t_win_val;
	win.start_t_win_val = stat.q[2] - win.t_bandwidth;
	win.end_t_win_val = stat.q[2] + win.t_bandwidth;

	//Case 1:t_{q_n} - t_q <= bandwidth
	if (stat.incr_t <= win.t_bandwidth)
	{
		cur_index = win.start_t_win_pos;
		for (int i = cur_index; i < stat.n; i++)
		{
			if (stat.sorted_featureVector[i][2] > win.start_t_win_val)
			{
				win.start_t_win_pos = i;
				isStart = true;
			}

			if (isStart == true)
				break;

			D_set.push_back(i);
		}

		cur_index = win.center_t_win_pos;
		for (int i = cur_index; i < stat.n; i++)
		{
			if (stat.sorted_featureVector[i][2] > stat.q[2])
			{
				win.center_t_win_pos = i;
				isCenter = true;
			}

			if (isCenter == true)
				break;

			C_set.push_back(i);
		}

		cur_index = win.end_t_win_pos + 1;
		for (int i = cur_index; i < stat.n; i++)
		{
			if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
			{
				win.end_t_win_pos = i - 1;
				isEnd = true;
			}

			if (isEnd == true)
				break;

			I_set.push_back(i);
		}

		//update sliding window (Equation 9)
		for (int i = 0; i < (int)D_set.size(); i++)
		{
			id = D_set[i];
			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);

			stat.sliding_window_L[0] -= kernel_s_value;
			stat.sliding_window_L[1] -= stat.sorted_featureVector[id][2] * kernel_s_value;
		}

		for (int i = 0; i < (int)C_set.size(); i++)
		{
			id = C_set[i];
			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);

			stat.sliding_window_L[0] += kernel_s_value;
			stat.sliding_window_L[1] += stat.sorted_featureVector[id][2] * kernel_s_value;
			stat.sliding_window_R[0] -= kernel_s_value;
			stat.sliding_window_R[1] -= stat.sorted_featureVector[id][2] * kernel_s_value;
		}

		for (int i = 0; i < (int)I_set.size(); i++)
		{
			id = I_set[i];
			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);

			stat.sliding_window_R[0] += kernel_s_value;
			stat.sliding_window_R[1] += stat.sorted_featureVector[id][2] * kernel_s_value;
		}
	}

	//Case 2:t_{q_n} - t_q > bandwidth && t_{q_n}-t_q <= 2 * bandwidth
	if (stat.incr_t > win.t_bandwidth && stat.incr_t < 2 * win.t_bandwidth)
	{
		cur_index = win.center_t_win_pos;
		for (int i = cur_index; i < stat.n; i++)
		{
			if (stat.sorted_featureVector[i][2] > win.start_t_win_val)
			{
				win.start_t_win_pos = i;
				isStart = true;
			}

			if (isStart == true)
				break;

			S_set.push_back(i);
		}

		cur_index = win.end_t_win_pos + 1;
		for (int i = cur_index; i < stat.n; i++)
		{
			if (stat.sorted_featureVector[i][2] > stat.q[2])
			{
				win.center_t_win_pos = i;
				isCenter = true;
			}

			if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
			{
				win.end_t_win_pos = i - 1;
				isEnd = true;
			}

			if (isCenter == false)
				A_set.push_back(i);

			if (isEnd == false)
				I_set.push_back(i);
			else
				break;
		}

		//update sliding window (Equation 9)
		stat.sliding_window_L[0] = stat.sliding_window_R[0];
		stat.sliding_window_L[1] = stat.sliding_window_R[1];
		stat.sliding_window_R[0] = 0;
		stat.sliding_window_R[1] = 0;

		for (int i = 0; i < (int)S_set.size(); i++)
		{
			id = S_set[i];
			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);

			stat.sliding_window_L[0] -= kernel_s_value;
			stat.sliding_window_L[1] -= stat.sorted_featureVector[id][2] * kernel_s_value;
		}

		for (int i = 0; i < (int)A_set.size(); i++)
		{
			id = A_set[i];
			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);

			stat.sliding_window_L[0] += kernel_s_value;
			stat.sliding_window_L[1] += stat.sorted_featureVector[id][2] * kernel_s_value;
			stat.sliding_window_R[0] -= kernel_s_value;
			stat.sliding_window_R[1] -= stat.sorted_featureVector[id][2] * kernel_s_value;
		}

		for (int i = 0; i < (int)I_set.size(); i++)
		{
			id = I_set[i];
			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);

			stat.sliding_window_R[0] += kernel_s_value;
			stat.sliding_window_R[1] += stat.sorted_featureVector[id][2] * kernel_s_value;
		}
	}

	//Case 3:t_{q_n} - t_q > 2 * bandwidth
	if (stat.incr_t > 2 * win.t_bandwidth)
	{
		cur_index = win.end_t_win_pos + 1;
		for (int i = cur_index; i < stat.n; i++)
		{
			if (isStart == false)
			{
				if (stat.sorted_featureVector[i][2] > win.start_t_win_val)
				{
					win.start_t_win_pos = i;
					isStart = true;
				}
			}

			if (isCenter == false)
			{
				if (stat.sorted_featureVector[i][2] > stat.q[2])
				{
					win.center_t_win_pos = i;
					isCenter = true;
				}
			}

			if (isEnd == false)
			{
				if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
				{
					win.end_t_win_pos = i - 1;
					isEnd = true;
				}
			}

			if (isStart == true && isEnd == false)
				I_set.push_back(i);

			if (isEnd == true)
				break;
		}

		//update sliding window (Equation 9)
		stat.sliding_window_L[0] = 0; stat.sliding_window_L[1] = 0;
		stat.sliding_window_R[0] = 0; stat.sliding_window_R[1] = 0;
		for (int i = 0; i < (int)I_set.size(); i++)
		{
			id = I_set[i];
			kernel_s_value = spatial_kernel(stat.q, stat.sorted_featureVector[id], stat);
			if (id < win.center_t_win_pos)
			{
				stat.sliding_window_L[0] += kernel_s_value;
				stat.sliding_window_L[1] += stat.sorted_featureVector[id][2] * kernel_s_value;
			}

			if (id >= win.center_t_win_pos)
			{
				stat.sliding_window_R[0] += kernel_s_value;
				stat.sliding_window_R[1] += stat.sorted_featureVector[id][2] * kernel_s_value;
			}
		}
	}

	win_density_value = (stat.sliding_window_L[0] + stat.sliding_window_R[0])
		- stat.gamma_t*(stat.q[2] * stat.sliding_window_L[0] - stat.sliding_window_L[1]
			+ stat.sliding_window_R[1] - stat.q[2] * stat.sliding_window_R[0]);

	return win_density_value;
}

void SWS(int x_index, int y_index, statistics& stat, double s_bandwidth, double t_bandwidth)
{
	win_status win;

	win.s_bandwidth = s_bandwidth; 
	win.t_bandwidth = t_bandwidth;

	obtain_q(x_index, y_index, 0, stat);
	if (stat.k_type_t == 0) //Triangular kernel
		stat.out_tensor[x_index][y_index][0] = compute_init_window_density_tri(stat, win);
	else //Epanechnikov and quartic kernels
		stat.out_tensor[x_index][y_index][0] = compute_init_window_density(stat, win);

	//Incremental update algorithm (Section 3.2 in our paper)
	for (int t_index = 1; t_index < stat.n_t; t_index++)
	{
		stat.q[2] = stat.t_L + t_index * stat.incr_t;

		if (stat.k_type_t == 0) //Triangular kernel
			stat.out_tensor[x_index][y_index][0] = incr_update_window_density_tri(stat, win);
		else //Epanechnikov and quartic kernels
			stat.out_tensor[x_index][y_index][t_index] = incr_update_window_density(stat, win);
	}
}