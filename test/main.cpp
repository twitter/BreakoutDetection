#include <vector>
#include <iostream>
#include "EdmPer.h"
#include "EdmMulti.h"
#include "EdmTail.h"
#include "Edmx.h"

using namespace std;
using namespace BreakoutDetection;

int main() {
    int min_size = 24;
    std::vector<double> Z;

    for (int j = 0; j < 10; j++) {
        for (int i = 0; i < min_size * 2; i++) {
            Z.push_back(0.0);
        }
        for (int i = 0; i < min_size * 2; i++) {
            Z.push_back(100.0);
        }
    }

    {
        double percent = 0;
        int degree = 0;

        EdmPer edmPer;
        edmPer.evaluate(Z, min_size, percent, degree);

        std::cout << "edmPer" << std::endl;
        for (int i = 0; i < edmPer.getLoc().size(); i++) {
            std::cout << edmPer.getLoc()[i] << std::endl;
        }
    }

    {
        double beta = 0;
        int degree = 0;

        EdmMulti edmMulti;
        edmMulti.evaluate(Z, min_size, beta, degree);

        std::cout << "edmMulti" << std::endl;
        for (int i = 0; i < edmMulti.getLoc().size(); i++) {
            std::cout << edmMulti.getLoc()[i] << std::endl;
        }
    }

    {
        double alpha = 2.0;
        double quant = 0.5;

        EdmTail edmTail;
        edmTail.evaluate(Z, min_size, alpha, quant);

        std::cout << "edmTail" << std::endl;
        std::cout << "loc: " << edmTail.getLoc() << std::endl;
        std::cout << "tail: " << edmTail.getTail() << std::endl;
        std::cout << "stat: " << edmTail.getStat() << std::endl;
    }

    {
        double alpha = 2.0;

        Edmx edmx;
        edmx.evaluate(Z, min_size, alpha);

        std::cout << "edmx" << std::endl;
        std::cout << "loc: " << edmx.getLoc() << std::endl;
        std::cout << "tail: " << edmx.getTail() << std::endl;
        std::cout << "stat: " << edmx.getStat() << std::endl;
    }

    return 0;
}