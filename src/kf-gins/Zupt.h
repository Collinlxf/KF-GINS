#pragma once

#include <array>
#include <vector>
#include <cmath>

#include "common/types.h"
#include "kf-gins/kf_gins_types.h"

class ZuptProcessor {
 public:
        explicit ZuptProcessor(int imudatarate);

        void AddData(const IMU& imuraw);
        bool IsZupt(const LcConfig& lcconfig, std::array<double, 6>* avg, double* quality);

 private:
    std::vector<double> zuptbuff_;
    int zuptindex_;
    int imudatarate_;
};
