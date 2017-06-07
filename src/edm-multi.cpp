#include <set>
#include <cmath>
#include <vector>
#include"helper.h"
#include "EdmMulti.h"

using namespace std;
using namespace BreakoutDetection;

//Z: time series
//min_size: minimum segment size
//beta: penalization term for the addition of a change point

void EdmMulti::evaluate(const std::vector<double> &Z, int min_size, double beta, int degree) {

	//identify which type of penalization to use
	double (*G)(double);
	switch (degree) {
		case 1:
			G = Linear;
			break;
		case 2:
			G = Quadratic;
			break;
		default:
			G = Const;
			break;
	}

	int n = Z.size();
	if (beta < 0)//assume that beta is a positive number
		beta = -beta;

	loc.clear();
	prev.clear();
	number.clear();
	F.clear();

	prev.resize(n + 1, 0);//store optimal location of previous change point
	number.resize(n + 1, 0);//store the number of change points in optimal segmentation
	F.resize(n + 1, -3);//store optimal statistic value
	//F[s] is calculated using observations { Z[0], Z[1], ..., Z[s-1] }

	//trees used to store the "upper half" of the considered observations
	multiset<double> right_min, left_min;
	//trees used to store the "lower half" of the considered observations
	multiset<double, std::greater<double> > right_max, left_max;

	//Iterate over possible locations for the last change
	for (int s = 2 * min_size; s < n + 1; ++s) {
		right_max.clear();
		right_min.clear();//clear right trees
		left_max.clear();
		left_min.clear();//clear left trees

		//initialize left and right trees to account for minimum segment size
		for (int i = prev[min_size - 1]; i < min_size - 1; ++i)
			insert_element(left_min, left_max, Z[i]);
		for (int i = min_size - 1; i < s; ++i)
			insert_element(right_min, right_max, Z[i]);

		//Iterate over possible locations for the penultiamte change
		for (int t = min_size; t < s - min_size + 1; ++t) {//modify limits to deal with min_size
			insert_element(left_min, left_max, Z[t - 1]);//insert element into left tree
			remove_element(right_min, right_max, Z[t - 1]);//remove element from right tree
			//left tree now has { Z[prev[t-1]], ..., Z[t-1] }
			//right tree now has { Z[t], ..., Z[s-1] }

			//check to see if optimal position of previous change point has changed
			//if so update the left tree
			if (prev[t] > prev[t - 1]) {
				for (int i = prev[t - 1]; i < prev[t]; ++i)
					remove_element(left_min, left_max, Z[i]);
			}
			else if (prev[t] < prev[t - 1]) {
				for (int i = prev[t]; i < prev[t - 1]; ++i)
					insert_element(left_min, left_max, Z[i]);
			}

			//calculate statistic value
			double left_median = get_median(left_min, left_max), right_median = get_median(right_min, right_max);
			double normalize = ((t - prev[t]) * (s - t)) / (std::pow((double) (s - prev[t]), 2.0));
			double tmp = F[t] + normalize * std::pow(left_median - right_median, 2.0) - beta * G(number[t]);
			//check for improved optimal statistic value
			if (tmp > F[s]) {
				number[s] = number[t] + 1;
				F[s] = tmp;
				prev[s] = t;
			}
		}
	}

	//obtain list of optimal change point estimates
	int at = n;
	while (at) {
		if (prev[at])//don't insert 0 as a change point estimate
			loc.push_back(prev[at]);
		at = prev[at];
	}
	sort(loc.begin(), loc.end());
}
