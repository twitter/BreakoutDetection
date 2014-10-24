#include<Rcpp.h>
#include<set>
#include<algorithm>
#include<cmath>

using namespace Rcpp;
using namespace std;

extern double Linear(double x){ return 1;}
extern double Const(double x){ return 0;}
extern double Quadratic(double x){ return 2*x+1;}


/*
Use 2 multisets (red-black trees) to keep track of the median. One tree for the larger (m) and 
one for the smaller (M) observations. Insertion and deletion in O(log(n)) and find 
the median in O(1), additional memory use is O(n).
*/

//insert x into the appropriate tree
extern void insert_element(multiset<double>& m, multiset<double, std::greater<double> >& M, double x){
	
	if(m.empty() || x < *(m.begin()))
		M.insert(x);
	else
		m.insert(x);
	if(m.size() > M.size() + 1){
		multiset<double>::iterator i;
		i = m.begin();
		m.erase(m.begin());
		M.insert(*i);
	}
	else if(M.size() > m.size() + 1){
		multiset<double, std::greater<double> >::iterator i;
		i = M.begin();
		M.erase(M.begin());
		m.insert(*i);
	}
}

//given a pair of trees obtain the median
extern double get_median(const multiset<double>& m, const multiset<double, std::greater<double> >& M){

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
extern void remove_element(multiset<double>& m, multiset<double, std::greater<double> >& M, double x){

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
		m.erase(m.begin());
		M.insert(*i);
	}
	else if(M.size() > m.size() + 1){
		multiset<double, std::greater<double> >::iterator i;
		i = M.begin();
		M.erase(M.begin());
		m.insert(*i);
	}
}
