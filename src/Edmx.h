#ifndef BREAKOUT_DETECTION_EDMX_H
#define BREAKOUT_DETECTION_EDMX_H

#include <vector>

namespace BreakoutDetection {

    class Edmx {
        public:
            void evaluate(const std::vector<double> &Z, int min_size = 24, double alpha = 2);

            double getLoc() {
                return loc;
            }
            double getTail() {
                return tail;
            }
            double getStat() {
                return stat;
            }

        private:
            double loc;
            double tail;
            double stat;
    };
}

#endif