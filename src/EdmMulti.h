#ifndef BREAKOUT_DETECTION_EDM_MULTI_H
#define BREAKOUT_DETECTION_EDM_MULTI_H

#include <vector>

namespace BreakoutDetection {

    class EdmMulti {
        public:
            void evaluate(const std::vector<double> &Z, int min_size = 24, double beta = 0, int degree = 0);

            std::vector<int> getLoc() {
                return loc;
            }
            std::vector<double> getF() {
                return F;
            }
            std::vector<int> getNumber() {
                return number;
            }
            std::vector<int> getPrev() {
                return prev;
            }

        private:
            std::vector<int> loc;
            std::vector<double> F;
            std::vector<int> number;
            std::vector<int> prev;
    };
}

#endif