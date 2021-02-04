#pragma once
#ifndef SWS_H
#define SWS_H

#include "SS.h"

struct index_time_pair
{
	int index;
	double time;
	bool operator<(const index_time_pair& pair) { return time < pair.time; }
};

struct win_status
{
	double s_bandwidth; 
	double t_bandwidth;

	//prev window
	double start_t_win_val_prev;
	double end_t_win_val_prev;

	//cur window
	double start_t_win_val; 
	double end_t_win_val;
	int start_t_win_pos;
	int	end_t_win_pos;
	int center_t_win_pos; //Used in triangular kernel
};

//void sort_vector_temporal(statistics& stat);
void init_SWS(statistics& stat);
double compute_init_window_density(statistics& stat, win_status& win);
void update_sliding_window(statistics& stat, vector<int>& index_set, bool is_positive);
double incr_update_window_density(statistics& stat, win_status& win);

double compute_init_window_density_tri(statistics& stat, win_status& win);
double incr_update_window_density_tri(statistics& stat, win_status& win);

void SWS(int x_index, int y_index, statistics& stat, double s_bandwidth, double t_bandwidth);

#endif