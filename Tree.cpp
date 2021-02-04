#include "Tree.h"

void Node::report_subtree(vector<int>& range_result_idList)
{
	for (int r = 0; r < (int)idList.size(); r++)
		range_result_idList.push_back(idList[r]);
}