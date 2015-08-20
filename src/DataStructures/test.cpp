#include <vector>
#include <iostream>

int main()
{
	using namespace std;

	vector<int> testVec(10, 1);

	testVec.assign(3, 2);

	for(auto it = testVec.begin(); it != testVec.end(); ++it)
	{
		cout << *it << endl;
	}

	return 0;
}