#ifndef BREAKOUT_DETECTION_HELPER_H
#define BREAKOUT_DETECTION_HELPER_H

#include<set>
#include<algorithm>
#include<cmath>

namespace BreakoutDetection {

    double get_median(std::multiset<double> &, std::multiset<double, std::greater<double> > &);
    void insert_element(std::multiset<double> &, std::multiset<double, std::greater<double> > &, double);
    void remove_element(std::multiset<double> &, std::multiset<double, std::greater<double> > &, double);

    double Linear(double x);
    double Const(double x);
    double Quadratic(double x);
}

#endif