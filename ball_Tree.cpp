#include "ball_Tree.h"

void ballNode::update_Augment(double**featureVector, int dim)
{
	int id;
	int num_of_points;
	double cur_radius_spatial;
	double cur_radius_temporal;

	center = new double[dim];
	for (int d = 0; d < dim; d++)
		center[d] = 0;

	num_of_points = (int)idList.size();
	//obtain the center for this node
	for (int i = 0; i < num_of_points; i++)
	{
		id = idList[i];
		for (int d = 0; d < dim; d++)
			center[d] += featureVector[id][d] / num_of_points;
	}

	//obtain the radius of this node
	radius_spatial = 0;
	radius_temporal = 0;
	for (int i = 0; i < num_of_points; i++)
	{
		id = idList[i];

		cur_radius_spatial = sqrt((center[0] - featureVector[id][0])*(center[0] - featureVector[id][0])
			+ (center[1] - featureVector[id][1])*(center[1] - featureVector[id][1]));
		cur_radius_temporal = fabs(center[2] - featureVector[id][2]);

		radius_spatial = max(radius_spatial, cur_radius_spatial);
		radius_temporal = max(radius_temporal, cur_radius_temporal);
	}
}

void ballNode::check_leaf(double*q, vector<int>& range_result_idList, double**featureVector, double s_bandwidth, double t_bandwidth)
{
	bool is_inside;
	int id;
	double spatial_dist, temporal_dist;

	for (int i = 0; i < (int)idList.size(); i++)
	{
		is_inside = true;
		id = idList[i];

		spatial_dist = s_dist(q, featureVector[id]);
		temporal_dist = t_dist(q, featureVector[id]);

		if (spatial_dist > s_bandwidth || temporal_dist > t_bandwidth)
			continue;

		range_result_idList.push_back(id);
	}
}

ballNode*ballNode::createNode()
{
	return new ballNode();
}

ballTree::ballTree(int dim, double**featureVector, int leafCapacity)
{
	this->dim = dim;
	this->featureVector = featureVector;
	this->leafCapacity = leafCapacity;
}

void ballTree::divide_node(ballNode*node, ballNode*leftNode, ballNode*rightNode)
{
	int index_1, index_2;
	int id, id_1, id_2;
	int best_id_1, best_id_2;
	double cur_dist;
	double max_dist = 0;
	int num_of_points;
	num_of_points = (int)node->idList.size();
	int half_size = (int)floor(num_of_points / 2.0);

	for (int r = 0; r < rand_num; r++)
	{
		index_1 = rand() % num_of_points; index_2 = rand() % num_of_points;
		id_1 = node->idList[index_1]; id_2 = node->idList[index_2];
		cur_dist = e_dist(featureVector[id_1], featureVector[id_2]);

		if (cur_dist > max_dist)
		{
			max_dist = cur_dist;
			best_id_1 = id_1;
			best_id_2 = id_2;
		}
	}

	for (int i = 0; i < num_of_points; i++)
	{
		id = node->idList[i];
		if ((int)leftNode->idList.size() >= half_size)
		{
			rightNode->idList.push_back(id);
			continue;
		}
		if ((int)rightNode->idList.size() >= half_size)
		{
			leftNode->idList.push_back(id);
			continue;
		}

		if (e_dist(featureVector[id], featureVector[id_1]) < e_dist(featureVector[id], featureVector[id_2]))
			leftNode->idList.push_back(id);
		else
			rightNode->idList.push_back(id);
	}
}

void ballTree::ballTree_Recur(ballNode*node)
{
	ballNode*leftNode;
	ballNode*rightNode;
	int num_of_points;

	num_of_points = (int)node->idList.size();
	//base case
	if (num_of_points <= leafCapacity)
		return;

	//create two children
	leftNode = node->createNode();
	rightNode = node->createNode();
	divide_node(node, leftNode, rightNode);
	leftNode->update_Augment(featureVector,dim);
	rightNode->update_Augment(featureVector, dim);

	ballTree_Recur(leftNode);
	ballTree_Recur(rightNode);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void ballTree::build_ballTree(statistics& stat)
{
	for (int i = 0; i < stat.n; i++)
		rootNode->idList.push_back(i);
	((ballNode*)rootNode)->update_Augment(featureVector, dim);

	ballTree_Recur((ballNode*)rootNode);
}

void ballTree::range_search_Recur(ballNode*node, double*q, double s_bandwidth, double t_bandwidth)
{
	double spatial_dist;
	double temporal_dist;
	ballNode*childNode;

	//base case
	if ((int)node->childVector.size() == 0)
	{
		node->check_leaf(q, range_result_idList, featureVector, s_bandwidth, t_bandwidth);
		return;
	}

	for (int c = 0; c < (int)node->childVector.size(); c++)
	{
		childNode = (ballNode*)node->childVector[c];
		spatial_dist = s_dist(q, childNode->center);
		temporal_dist = t_dist(q, childNode->center);

		if (spatial_dist + childNode->radius_spatial <= s_bandwidth && temporal_dist + childNode->radius_temporal <= t_bandwidth)
		{
			childNode->report_subtree(range_result_idList);
			return;
		}
			
		if (spatial_dist >= s_bandwidth + childNode->radius_spatial || temporal_dist >= t_bandwidth + childNode->radius_temporal)
			continue;
		
		range_search_Recur(childNode, q, s_bandwidth, t_bandwidth);
	}
}

void ballTree::range_search(double*q, double s_bandwidth, double t_bandwidth)
{
	double spatial_dist, temporal_dist;

	spatial_dist = s_dist(q, ((ballNode*)rootNode)->center);
	temporal_dist = t_dist(q, ((ballNode*)rootNode)->center);

	if (spatial_dist + ((ballNode*)rootNode)->radius_spatial <= s_bandwidth && temporal_dist + ((ballNode*)rootNode)->radius_temporal <= t_bandwidth)
		((ballNode*)rootNode)->report_subtree(range_result_idList);

	if (spatial_dist >= s_bandwidth + ((ballNode*)rootNode)->radius_spatial || temporal_dist >= t_bandwidth + ((ballNode*)rootNode)->radius_temporal)
		return;
	else
		range_search_Recur(((ballNode*)rootNode), q, s_bandwidth, t_bandwidth);
}