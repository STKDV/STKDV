#pragma once
#ifndef BALL_TREE_H
#define BALL_TREE_H

#include "Tree.h"

class ballNode : public Node
{
public:
	double*center;
	double radius_spatial;
	double radius_temporal;

	void update_Augment(double**featureVector, int dim);
	void check_leaf(double*q, vector<int>& range_result_idList, double**featureVector, double s_bandwidth, double t_bandwidth);
	ballNode*createNode();
};

class ballTree : public Tree
{
public:
	const int rand_num = 5;

	//offline
	ballTree(int dim, double**featureVector, int leafCapacity);
	void divide_node(ballNode*node, ballNode*leftNode, ballNode*rightNode);
	void ballTree_Recur(ballNode*node);
	void build_ballTree(statistics& stat);

	//online
	void range_search_Recur(ballNode*node, double*q, double s_bandwidth, double t_bandwidth);
	void range_search(double*q, double s_bandwidth, double t_bandwidth);
};

#endif