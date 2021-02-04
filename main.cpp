#include "alg_visual.h"

int main(int argc, char**argv)
{
	statistics stat;
	init_visual(argc, argv, stat);
	visual_Algorithm(stat);
	save_result_to_file(stat);
}