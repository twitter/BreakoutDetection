#include <vector>
#include <Rcpp.h>
#include "EdmPer.h"
#include "EdmMulti.h"
#include "EdmTail.h"
#include "Edmx.h"

using namespace Rcpp;

std::vector<double> toNumericVector(const Rcpp::NumericVector& Z) {
    std::vector<double> z;

    for (Rcpp::NumericVector::const_iterator it = Z.begin(); it != Z.end(); ++it) {
        z.push_back(*it);
    }
    return z;
}

// [[Rcpp::export]]
Rcpp::List EDM_percent(const Rcpp::NumericVector& Z, int min_size = 24, double percent = 0, int degree = 0) {

    BreakoutDetection::EdmPer edmPer;
    std::vector<double> z = toNumericVector(Z);
    edmPer.evaluate(z, min_size, percent, degree);

    //return statment for debugging
    //return List::create(_["loc"]=edmPer.getLoc(), _["F"]=edmPer.getF(), _["number"]=edmPer.getNumber,_["prev"]=edmPer.getPrev());

    return Rcpp::List::create(_["loc"]=edmPer.getLoc());
}

// [[Rcpp::export]]
Rcpp::List EDM_multi(const Rcpp::NumericVector& Z, int min_size = 24, double beta = 0, int degree = 0) {

    BreakoutDetection::EdmMulti edmMulti;
    std::vector<double> z = toNumericVector(Z);
    edmMulti.evaluate(z, min_size, beta, degree);

    //return statment for debugging
    //return List::create(_["loc"]=edmMulti.getLoc(), _["F"]=edmMulti.getF(), _["number"]=edmMulti.getNumber,_["prev"]=edmMulti.getPrev());

    return Rcpp::List::create(_["loc"]=edmMulti.getLoc());
}

// [[Rcpp::export]]
Rcpp::List EDM_tail(const Rcpp::NumericVector& Z, int min_size, double alpha, double quant) {

    BreakoutDetection::EdmTail edmTail;
    std::vector<double> z = toNumericVector(Z);
    edmTail.evaluate(z, min_size, alpha, quant);

    //return statment for debugging
    //return List::create(_["loc"]=edmTail.getLoc(), _["tail"]=edmTail.getTail(), _["stat"]=edm.getStat());

    return Rcpp::List::create(_["loc"]=edmTail.getLoc(), _["stat"]=edmTail.getStat());
}

// [[Rcpp::export]]
Rcpp::List EDMX(const Rcpp::NumericVector& Z, int min_size, double alpha) {

    BreakoutDetection::Edmx edmx;
    std::vector<double> z = toNumericVector(Z);
    edmx.evaluate(z, min_size, alpha);

    //return statment for debugging
    //return List::create(_["loc"]=edmx.getLoc(), _["tail"]=edmx.getTail(), _["stat"]=edmx.getStat());

    return Rcpp::List::create(_["loc"]=edmx.getLoc(), _["stat"]=edmx.getStat());
}
