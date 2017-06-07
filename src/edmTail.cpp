/*
This version calculates the between distance using the delta points around the change point estimate.
*/

#include<cmath>
#include<algorithm>
#include<iostream>
#include<vector>
#include "EdmTail.h"
#include"helper.h"

using namespace std;
using namespace BreakoutDetection;


//Class used to hold all the information about the
//breakout location and the interval trees
struct Information{
	std::vector<double> A, B, AB;
	double best_stat;
	int best_loc, best_t2;
	int min_size, b;
	
	Information(int, int);
};

Information::Information(int bb, int m){
	A = std::vector<double>(1<<(bb+1));
	B = std::vector<double>(1<<(bb+1));
	AB = std::vector<double>(1<<(bb+1));
	b = bb;
	best_stat = best_loc = best_t2 = -3;
	min_size = m;
}

/*void printInformation(Information& info){
	std::cout<<"best_stat: "<<info.best_stat<<std::endl;
	std::cout<<"best_t2: "<<info.best_t2<<std::endl;
	std::cout<<"best_loc: "<<info.best_loc<<std::endl;
}*/

void BackwardUpdate(const std::vector<double>& Z, Information& info, int& tau1, double quant, double alpha);
void ForwardUpdate(const std::vector<double>& Z, Information& info, int& tau1, double quant, double alpha);

int GetIndex(int B, double x){
	//Get index of leaf node interval containing x
	return (int)std::ceil(std::abs(x) * (1<<B)) + (1<<B) - 1;
}

double GetQuantile(std::vector<double>& x, double quant){
	//Return approximate quantile based on the interval tree

	int N = x.size();
	int k = std::ceil(x[1]*quant);
	double l=0, u=1;
	int i=1,j;
	while(i < N){ //Make sure that we do not go beyond the array bounds
		j = i<<1;
		if(j >= N)
			break;
		if(x[i] == k){//Exactly k elements in this node's subtree. So can terminate early
			//Return a weighted combination of the child node medians
			double lWeight = x[j]/(x[j]+x[j+1]);
			double rWeight = 1 - lWeight;
			double lu, rl;
			lu = (u+l)/2;
			rl = (u+lu)/2;
			return lWeight*(quant*(lu-l)+l) + rWeight*(quant*(u-rl)+rl);
		}
		else if(x[j] >= k){//More than k elements in node's left child's subtree, move to left child
			i = j;
			u = (l+u)/2;
		}
		else if(x[j] < k){//Not enough elements in node's left child's subtree, move to right child
			k -= x[j];
			i = j+1;
			l = (l+u)/2;
		}
	}
	return quant*(u-l)+l;
}

std::vector<int> AddToTree(int B, std::vector<double>& x){
	std::vector<int> A(1<<(B+1));
	std::vector<double>::iterator i;
	for(i = x.begin(); i < x.end(); ++i){//Iterage over items we wish to add to the tree
		int index = GetIndex(B,*i);
		while(index){
			++A[index];
			index /= 2;
		}
	}
	return A;
}

void EdmTail::evaluate(const std::vector<double>& Z, int min_size, double alpha, double quant){

	int N = Z.size();
	int eps = (int)std::ceil( std::log(N) );
	eps = std::max( eps, 10 );
	
	Information info(eps,min_size);

	int tau1 = info.min_size;
	int tau2 = tau1 * 2;

	
	
	//Populate trees and calculate statistic value for starting configuration of 
	//2 min_size segments
	for(int i=0; i<tau1; ++i){
		for(int j=i+1; j<tau1; ++j){
			int index = GetIndex(info.b,Z[i]-Z[j]);
			while(index){
				++info.A[index];
				index /= 2;
			}
		}
	}

	//Populate trees and calculate statistic value for starting configuration of 
	//2 min_size segments
	for(int i=tau1; i<tau2; ++i){
		for(int j=i+1; j<tau2; ++j){
			int index = GetIndex(info.b,Z[i]-Z[j]);
			while(index){
				++info.B[index];
				index /= 2;
			}
		}
	}

	//Populate trees and calculate statistic value for starting configuration of 
	//2 min_size segments
	for(int i=0; i<tau1; ++i){
		for(int j=tau1; j<tau2; ++j){
			int index = GetIndex(info.b,Z[i]-Z[j]);
			while(index){
				++info.AB[index];
				index /= 2;
			}
		}
	}

	double qa, qb, qc, stat;

	qa = std::pow( GetQuantile(info.A,quant), alpha);
	qb = std::pow( GetQuantile(info.B,quant), alpha);
	qc = std::pow( GetQuantile(info.AB,quant), alpha);

	stat = 2*qc - qa - qb;
	stat *= (double)(tau1)*(tau2-tau1)/(tau2);

	info.best_stat = stat;
	info.best_loc = tau1;
	info.best_t2 = tau2;

	//Increment tau2 and update trees and statistic 
	++tau2;
	for(; tau2<N+1; ++tau2){
		int index = GetIndex(info.b,Z[tau2-1]-Z[tau2-2]);
		while(index){//array position 0 is not used, so we exit once we reach this location
			++info.B[index];
			index /= 2;
		}
		qb = std::pow( GetQuantile(info.B,quant), alpha);
		stat = 2*qc - qa - qb;
		stat *= (double)(tau2-tau1)*tau1/tau2;

		if(stat > info.best_stat){
			info.best_stat = stat;
			info.best_loc = tau1;
			info.best_t2 = tau2;
		}
	}

	bool forward_move = false;
	//Initial consideration of other possible locations for tau1
	while(tau1 < N-min_size){
		//"warm start" to update tree and statistic value for other prefix series
		if(forward_move){
			ForwardUpdate(Z, info, tau1, quant, alpha);
		}
		else{
			BackwardUpdate(Z, info, tau1, quant, alpha);
		}
		forward_move = !forward_move;
	}

	this->loc = info.best_loc;
	this->tail = info.best_t2;
	this->stat = info.best_stat;
}

void ForwardUpdate(const std::vector<double>& Z, Information& info, int& tau1, double quant, double alpha){
	
	int min_size = info.min_size;
	int tau2 = tau1 + min_size;
	++tau1;
	int N = Z.size(), index;
	//Update A tree
	for(int i=tau1-min_size; i<tau1-1; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-1]);
		while(index){
			++info.A[index];
			index /= 2;
		}
	}
	for(int i=tau1-min_size; i<tau1; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-min_size-1]);
		while(index){
			--info.A[index];
			index /= 2;
		}
	}
	index = GetIndex(info.b, Z[tau1-min_size-1]-Z[tau1-min_size]);
	while(index){
		++info.A[index];
		index /= 2;
	}
	double qa = std::pow( GetQuantile(info.A,quant), alpha );

	//Update AB tree
	index = GetIndex(info.b, Z[tau1-1]-Z[tau1-min_size-1]);
	while(index){
		--info.AB[index];
		index /= 2;
	}
	for(int i=tau1; i<tau2; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-min_size-1]);
		while(index){
			--info.AB[index];
			index /= 2;
		}
		index = GetIndex(info.b, Z[i]-Z[tau1-1]);
		while(index){
			++info.AB[index];
			index /= 2;
		}
	}
	for(int i=tau1-min_size; i<tau1-1; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-1]);
		while(index){
			--info.AB[index];
			index /= 2;
		}
		index = GetIndex(info.b, Z[i]-Z[tau2]);
		while(index){
			++info.AB[index];
			index /= 2;
		}
	}
	index = GetIndex(info.b, Z[tau1-1]-Z[tau2]);
	while(index){
		++info.AB[index];
		index /= 2;
	}
	double qc = std::pow( GetQuantile(info.AB,quant), alpha );

	//Update B tree
	for(int i=tau1; i<tau2; ++i){
		index = GetIndex(info.b,Z[i]-Z[tau1-1]);
		while(index){
			--info.B[index];
			index /= 2;
		}
		index = GetIndex(info.b,Z[i]-Z[tau2]);
		while(index){
			++info.B[index];
			index /= 2;
		}
	}

	//Increment tau2 and update statistic value as we proceed
	++tau2;
	for(; tau2<N+1; ++tau2){
		index = GetIndex(info.b, Z[tau2-1]-Z[tau2-2]);
		while(index){
			++info.B[index];
			index /= 2;
		}
		double qb = std::pow( GetQuantile(info.B,quant), alpha );

		double stat = 2*qc - qa - qb;
		stat *= (double)(tau2-tau1)*tau1/tau2;

		if(stat > info.best_stat){
			info.best_stat = stat;
			info.best_loc = tau1;
			info.best_t2 = tau2;
		}
	}
}

void BackwardUpdate(const std::vector<double>& Z, Information& info, int& tau1, double quant, double alpha){

	int min_size = info.min_size;
	int tau2 = tau1 + min_size;
	++tau1;
	int N = Z.size(), index;
	//Update A tree
	for(int i=tau1-min_size; i<tau1-1; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-1]);
		while(index){
			++info.A[index];
			index /= 2;
		}
	}
	for(int i=tau1-min_size; i<tau1; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-min_size-1]);
		while(index){
			--info.A[index];
			index /= 2;
		}
	}
	index = GetIndex(info.b, Z[tau1-min_size-1]-Z[tau1-min_size]);
	while(index){
		++info.A[index];
		index /= 2;
	}
	double qa = std::pow( GetQuantile(info.A,quant), alpha );

	//Update AB tree
	index = GetIndex(info.b, Z[tau1-1]-Z[tau1-min_size-1]);
	while(index){
		--info.AB[index];
		index /= 2;
	}
	for(int i=tau1; i<tau2; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-min_size-1]);
		while(index){
			--info.AB[index];
			index /= 2;
		}
		index = GetIndex(info.b, Z[i]-Z[tau1-1]);
		while(index){
			++info.AB[index];
			index /= 2;
		}
	}
	for(int i=tau1-min_size; i<tau1-1; ++i){
		index = GetIndex(info.b, Z[i]-Z[tau1-1]);
		while(index){
			--info.AB[index];
			index /= 2;
		}
		index = GetIndex(info.b, Z[i]-Z[tau2]);
		while(index){
			++info.AB[index];
			index /= 2;
		}
	}
	index = GetIndex(info.b, Z[tau1-1]-Z[tau2]);
	while(index){
		++info.AB[index];
		index /= 2;
	}
	double qc = std::pow( GetQuantile(info.AB,quant), alpha );


	//Update B tree
	for(int i=tau1; i<tau1+min_size-1; ++i){
		index = GetIndex(info.b, Z[tau1+min_size-1]-Z[i]);
		while(index){
			++info.B[index];
			index /= 2;
		}
		index = GetIndex(info.b, Z[i]-Z[tau1-1]);
		while(index){
			--info.B[index];
			index /= 2;
		}
	}
	double qb = std::pow( GetQuantile(info.B,quant), alpha);
	//Move tau2 from the end of the time series to the front.
	//Update the statistic value along the way
	tau2 = N;

	for(; tau2>=tau1+min_size; --tau2){
		index = GetIndex(info.b,Z[tau2-1]-Z[tau2-2]);
		while(index){
			--info.B[index];
			index /= 2;
		}
		qb = std::pow( GetQuantile(info.B,quant), alpha);

		double stat = 2*qc - qa - qb;
		stat *= (double)(tau2-tau1)*tau1/tau2;

		if(stat > info.best_stat){
			info.best_stat = stat; 
			info.best_loc = tau1;
			info.best_t2 = tau2;
		}
	}
}
