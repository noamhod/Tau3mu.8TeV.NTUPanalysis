#ifndef RUNINFO_H
#define RUNINFO_H

#include "std.h"
// #include <algorithm>

//// Can use class, struct or primitive data type for Item=RunInfo
struct RunInfo
{
	int evt;
	int lbn;
	int run;
	
	bool operator==(const RunInfo& m) const
	{
		return((m.evt==evt)&&(m.lbn==lbn)&&(m.run==run));
	}
};
typedef vector<RunInfo>         RunInfoVector;
typedef RunInfoVector::iterator RunInfoIterator;

//// make a global instance of the vector
RunInfoVector vRunInfo;

//// the vector needs to be cleared per run
void clearRunInfo()
{
	vRunInfo.clear();
}

bool isDuplicatedEvent(int evt, int lbn, int run)
{
	//// init data for RunInfoVector
	RunInfo thisRunInfo;
	thisRunInfo.evt = evt;
	thisRunInfo.lbn = lbn;
	thisRunInfo.run = run;
	
	//// try to find thisRunInfo
	RunInfoIterator itRunInfo = find(vRunInfo.begin(), vRunInfo.end(), thisRunInfo);
	
	//// thisRunInfo found
	if(itRunInfo != vRunInfo.end())
	{
		_INF(1,"Skipping a duplicated event: evt:lbn:run="<<evt<<":"<<lbn<<":"<<run);
		return true;
	}

	//// thisRunInfo is not found
	vRunInfo.push_back(thisRunInfo);
	return false;
}


/*------------------------------------------
//// http://stackoverflow.com/questions/571394/how-to-find-an-item-in-a-stdvector
http://stackoverflow.com/questions/11340802/find-struct-in-vector
#include <algorithm>
#include <vector>

// You can use class, struct or primitive data type for Item
struct Item {
    //Some fields
	// and the operator = definition
};
typedef std::vector<Item> ItemVector;
typedef ItemVector::iterator ItemIterator;
//...
ItemVector vtItem;
//... (init data for vtItem)
Item itemToFind;
//...

ItemIterator itemItr;
itemItr = std::find(vtItem.begin(), vtItem.end(), itemToFind);
if (itemItr != vtItem.end()) {
    // Item found
    // doThis()
}
else {
    // Item not found
    // doThat()
}
------------------------------------------*/

#endif