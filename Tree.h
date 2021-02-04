#pragma once
#ifndef TREE_H
#define TREE_H

#include "init_visual.h"

class Node;
class Tree;
const int leafCapacity = 40;

class Tree
{
public:
	int dim;
	double**featureVector;
	int leafCapacity; //Set this value to 40 in tKDC/Scikit-learn software
	Node*rootNode;
	vector<int> range_result_idList; //Used in range query

	//virtual void range_search(double**query_boundary) = 0;
};

class Node
{
public:
	vector<int> idList; //idxs
	vector<Node*> childVector;

	void report_subtree(vector<int>& range_result_idList);
	virtual Node*createNode() = 0;
	//virtual void update_Aug(Node*node, Tree*t) = 0;

	//Tree tree;
};

#endif