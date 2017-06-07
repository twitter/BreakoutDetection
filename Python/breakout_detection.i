%module breakout_detection

%{
#include "../src/EdmPer.h"
#include "../src/EdmMulti.h"
#include "../src/EdmTail.h"
#include "../src/Edmx.h"
%}

%include "std_vector.i"

namespace std {
    %template(IntegerVector) vector<int>;
    %template(NumericVector) vector<double>;
};

%include "../src/EdmPer.h"
%include "../src/EdmMulti.h"
%include "../src/EdmTail.h"
%include "../src/Edmx.h"