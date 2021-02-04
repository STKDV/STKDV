#include "alg_visual.h"

void visual_Algorithm(statistics& stat)
{
	kdTree kd_tree(stat.dim, stat.featureVector, leafCapacity);
	ballTree ball_tree(stat.dim, stat.featureVector, leafCapacity);
	double s_bandwidth, t_bandwidth;
	double run_time;

	auto start_s = chrono::high_resolution_clock::now();
	if (stat.method == 1) //RQS (kd-tree)
	{
		stat.query_boundary = new double*[stat.dim];
		for (int d = 0; d < stat.dim; d++)
			stat.query_boundary[d] = new double[2];
		kd_tree.rootNode = new kdNode();
		kd_tree.build_kdTree(stat);
	}
	
	if (stat.method == 2 || stat.method == 4) //SWS for method = 4, RQS (ball-tree) for method = 2
	{
		s_bandwidth = 1.0 / stat.gamma_s;
		t_bandwidth = 1.0 / stat.gamma_t;

		if (stat.method == 2) 
		{
			ball_tree.rootNode = new ballNode();
			ball_tree.build_ballTree(stat);
		}
		if (stat.method == 4)
			init_SWS(stat);
	}

	for (int x_index = 0; x_index < stat.n_x; x_index++)
	{
		for (int y_index = 0; y_index < stat.n_y; y_index++)
		{
			if (stat.method == 4) //SWS
			{
				SWS(x_index, y_index, stat, s_bandwidth, t_bandwidth);
				continue;
			}

			for (int t_index = 0; t_index < stat.n_t; t_index++)
			{
				obtain_q(x_index, y_index, t_index, stat);
				if (stat.method == 0) //SCAN method
					stat.out_tensor[x_index][y_index][t_index] = SCAN(stat);
				if (stat.method == 1) //RAQ with kd-tree
				{
					kd_tree.obtain_boundary(stat);
					kd_tree.range_search(stat.query_boundary);
					stat.out_tensor[x_index][y_index][t_index] = SCAN_RQS_result(stat, kd_tree);
					kd_tree.range_result_idList.clear();
				}
				if (stat.method == 2) //RAQ with ball-tree
				{
					ball_tree.range_search(stat.q, s_bandwidth, t_bandwidth);
					stat.out_tensor[x_index][y_index][t_index] = SCAN_RQS_result(stat, ball_tree);
					ball_tree.range_result_idList.clear();
				}
			}
		}
	}

	auto end_s = chrono::high_resolution_clock::now();
	run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;
	std::cout << "method " << stat.method << ":" << run_time << endl;

	save_result_to_file(stat);
}