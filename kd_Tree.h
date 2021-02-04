#pragma once
#ifndef KD_TREE_H
#define KD_TREE_H

#include "Tree.h"

class kdNode : public Node  //Gan_SIGMOD17
{
public:
	double**boundary; //boundary

	void update_boundary(Tree*t);
	int check_condition_boundaries(double**query_boundary, int dim); //Three conditions: Cover (0), Intersect (1) and Not intersect (2)
	//void report_subtree(vector<int>& range_result_idList);
	void check_leaf(double**query_boundary, double**featureVector, vector<int>& range_result_idList, int dim);
	kdNode*createNode();
};

class kdTree : public Tree
{
public:
	//offline phase
	kdTree(int dim, double**featureVector, int leafCapacity);
	double obtain_SplitValue(kdNode*node, int split_Dim);
	void KD_Tree_Recur(kdNode*node, int split_Dim);
	void build_kdTree(statistics& stat);

	//online phase
	void range_search_Recur(kdNode*node, double**query_boundary);
	void range_search(double**query_boundary);
	void obtain_boundary(statistics& stat);
};

#endif