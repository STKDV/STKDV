#include "kd_Tree.h"

void kdNode::update_boundary(Tree*t)
{
	double x_min, x_max;
	double y_min, y_max;
	double t_min, t_max;
	int id;

	x_min = inf; y_min = inf; t_min = inf;
	x_max = -inf; y_max = -inf; t_max = -inf;
	boundary = new double*[t->dim];
	for (int d = 0; d < t->dim; d++)
		boundary[d] = new double[2];
	
	for (int i = 0; i < (int)this->idList.size(); i++)
	{
		id = this->idList[i];

		if (t->featureVector[id][0] < x_min)
			x_min = t->featureVector[id][0];
		if (t->featureVector[id][0] > x_max)
			x_max = t->featureVector[id][0];
		if (t->featureVector[id][1] < y_min)
			y_min = t->featureVector[id][1];
		if (t->featureVector[id][1] > y_max)
			y_max = t->featureVector[id][1];
		if (t->featureVector[id][2] < t_min)
			t_min = t->featureVector[id][2];
		if (t->featureVector[id][2] > t_max)
			t_max = t->featureVector[id][2];
	}

	boundary[0][0] = x_min; boundary[0][1] = x_max;
	boundary[1][0] = y_min; boundary[1][1] = y_max;
	boundary[2][0] = t_min; boundary[2][1] = t_max;
}

//Three conditions: Cover (0), Intersect (1) and Not intersect (2)
int kdNode::check_condition_boundaries(double**query_boundary, int dim)
{
	//initialize to the value 0 (Cover)
	int condition = 0;
	int cur_condition = 2; //This state is "not intersect"

	for (int d = 0; d < dim; d++)
	{
		if (query_boundary[d][0] <= boundary[d][0] && boundary[d][1] <= query_boundary[d][1])
			cur_condition = 0;
		if (boundary[d][0] <= query_boundary[d][0] && query_boundary[d][1] <= boundary[d][1])
			cur_condition = 1;
		if (query_boundary[d][0] <= boundary[d][0] && boundary[d][0] <= query_boundary[d][1] && query_boundary[d][1] <= boundary[d][1])
			cur_condition = 1;
		if (boundary[d][0] <= query_boundary[d][0] && query_boundary[d][0] <= boundary[d][1] && boundary[d][1] <= query_boundary[d][1])
			cur_condition = 1;

		if (cur_condition == 2) //N + C => N, N + I => N, //N + N => N 
		{
			condition = 2;
			break;
		}

		if (cur_condition == 1 && (condition == 0 || condition == 1))
			condition = 1;
	}

	return condition;
}

/*void kdNode::report_subtree(vector<int>& range_result_idList)
{
	for (int r = 0; r < (int)idList.size(); r++)
		range_result_idList.push_back(idList[r]);
}*/

void kdNode::check_leaf(double**query_boundary, double**featureVector, vector<int>& range_result_idList, int dim)
{
	bool is_inside;
	int id;
	for (int i = 0; i < (int)idList.size(); i++)
	{
		is_inside = true;
		id = idList[i];
		for (int d = 0; d < dim; d++)
		{
			if (query_boundary[d][0] <= featureVector[id][d] && featureVector[id][d] <= query_boundary[d][1])
				continue;
			is_inside = false;
			break;
		}

		if (is_inside == true)
			range_result_idList.push_back(id);
	}
}

kdNode*kdNode::createNode()
{
	return new kdNode();
}

kdTree::kdTree(int dim, double**featureVector, int leafCapacity)
{
	this->dim = dim;
	this->featureVector = featureVector;
	this->leafCapacity = leafCapacity;
}

double kdTree::obtain_SplitValue(kdNode*node, int split_Dim)
{
	vector<double> tempVector;
	int id;
	int middle_left, middle_right, middle;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		tempVector.push_back(featureVector[id][split_Dim]);
	}

	sort(tempVector.begin(), tempVector.end());

	if ((int)tempVector.size() % 2 == 0)//even number
	{
		middle_right = (int)tempVector.size() / 2;
		middle_left = middle_right - 1;

		return ((tempVector[middle_left] + tempVector[middle_right]) / 2.0);
	}
	else
	{
		middle = ((int)tempVector.size() - 1) / 2;
		return tempVector[middle];
	}

	tempVector.clear();
}

void kdTree::KD_Tree_Recur(kdNode*node, int split_Dim)
{
	int id;
	int counter;
	kdNode*leftNode;
	kdNode*rightNode;
	double splitValue;

	//base case
	if ((int)node->idList.size() <= leafCapacity)
		return;

	splitValue = obtain_SplitValue(node, split_Dim); //code here

	//create two children
	leftNode = node->createNode();
	rightNode = node->createNode();

	counter = 0;
	int halfSize = ((int)node->idList.size()) / 2;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		if (featureVector[id][split_Dim] <= splitValue && counter <= halfSize)
		{
			leftNode->idList.push_back(id);
			counter++;
		}
		else
			rightNode->idList.push_back(id);
	}

	leftNode->update_boundary(this);
	rightNode->update_boundary(this);

	KD_Tree_Recur(leftNode, (split_Dim + 1) % dim);
	KD_Tree_Recur(rightNode, (split_Dim + 1) % dim);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void kdTree::build_kdTree(statistics& stat)
{
	for (int i = 0; i < stat.n; i++)
		rootNode->idList.push_back(i);
	((kdNode*)rootNode)->update_boundary(this);

	KD_Tree_Recur((kdNode*)rootNode, 0);
}

//online phase
//Follows the algorithm in the book (Computational Geometry: Algorithms and Applications (Second Edition) p.103)
void kdTree::range_search_Recur(kdNode*node, double**query_boundary)
{
	int condition;
	kdNode*child_node;

	if (node->childVector.size() == 0) //child node
		node->check_leaf(query_boundary, featureVector, range_result_idList, dim);
	else
	{
		for (int c = 0; c < (int)node->childVector.size(); c++)
		{
			child_node = (kdNode*)node->childVector[c];

			condition = child_node->check_condition_boundaries(query_boundary, dim);
			if (condition == 0) //Cover
				child_node->report_subtree(range_result_idList);
			else if (condition == 1)//Intersect
				range_search_Recur(child_node, query_boundary);
		}
	}
}

void kdTree::range_search(double**query_boundary)
{
	int condition;

	condition = ((kdNode*)rootNode)->check_condition_boundaries(query_boundary, dim);
	if (condition == 0)
		((kdNode*)rootNode)->report_subtree(range_result_idList);
	if (condition == 1)
		range_search_Recur((kdNode*)rootNode, query_boundary);
}

void kdTree::obtain_boundary(statistics& stat)
{
	stat.query_boundary[0][0] = stat.q[0] - (1.0 / stat.gamma_s);
	stat.query_boundary[0][1] = stat.q[0] + (1.0 / stat.gamma_s);
	stat.query_boundary[1][0] = stat.q[1] - (1.0 / stat.gamma_s);
	stat.query_boundary[1][1] = stat.q[1] + (1.0 / stat.gamma_s);
	stat.query_boundary[2][0] = stat.q[2] - (1.0 / stat.gamma_t);
	stat.query_boundary[2][1] = stat.q[2] + (1.0 / stat.gamma_t);
}