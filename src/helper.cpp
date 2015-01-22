#include<set>
#include<cmath>
#include "helper.h"

using namespace std;
using namespace BreakoutDetection;

double BreakoutDetection::Linear(double x){ return 1;}
double BreakoutDetection::Const(double x){ return 0;}
double BreakoutDetection::Quadratic(double x){ return 2*x+1;}


/*
Use 2 multisets (red-black trees) to keep track of the median. One tree for the larger (m) and 
one for the smaller (M) observations. Insertion and deletion in O(log(n)) and find 
the median in O(1), additional memory use is O(n).
*/

//insert x into the appropriate tree
void BreakoutDetection::insert_element(multiset<double>& m, multiset<double, std::greater<double> >& M, double x){
	
	if(m.empty() || x < *(m.begin()))
		M.insert(x);
	else
		m.insert(x);
	if(m.size() > M.size() + 1){
		multiset<double>::iterator i;
		i = m.begin();
		M.insert(*i);
		m.erase(m.begin());
	}
	else if(M.size() > m.size() + 1){
		multiset<double, std::greater<double> >::iterator i;
		i = M.begin();
		m.insert(*i);
		M.erase(M.begin());
	}
}

//given a pair of trees obtain the median
double BreakoutDetection::get_median(multiset<double>& m, multiset<double, std::greater<double> >& M){

	if(m.size() > M.size())
		return *(m.begin());
	else if(M.size() > m.size())
		return *(M.begin());
	else
		return ( *(M.begin()) + *(m.begin()) )/2;
}

//remove x from the tree, if multiple copies of x exist only remove 1
//since this method is never called by the user directly it is assumed 
//that there is at least 1 copy of x
void BreakoutDetection::remove_element(multiset<double>& m, multiset<double, std::greater<double> >& M, double x){

	if(x < *(m.begin())){
		multiset<double, std::greater<double> >::iterator i = M.find(x);
		M.erase(i);
	}
	else{
		multiset<double>::iterator i = m.find(x);
		m.erase(i);
	}
	if(m.size() > M.size() + 1){
		multiset<double>::iterator i;
		i = m.begin();
		M.insert(*i);
		m.erase(m.begin());
	}
	else if(M.size() > m.size() + 1){
		multiset<double, std::greater<double> >::iterator i;
		i = M.begin();
		m.insert(*i);
		M.erase(M.begin());
	}
}
