#include "Zupt.h"

ZuptProcessor::ZuptProcessor(int imudatarate)
    : zuptbuff_(imudatarate * 6, 0.0), zuptindex_(0), imudatarate_(imudatarate) {}

void ZuptProcessor::AddData(const IMU& imuraw) {
    zuptbuff_[zuptindex_ * 6 + 0] = imuraw.dtheta[0];
    zuptbuff_[zuptindex_ * 6 + 1] = imuraw.dtheta[1];
    zuptbuff_[zuptindex_ * 6 + 2] = imuraw.dtheta[2];
    zuptbuff_[zuptindex_ * 6 + 3] = imuraw.dvel[0];
    zuptbuff_[zuptindex_ * 6 + 4] = imuraw.dvel[1];
    zuptbuff_[zuptindex_ * 6 + 5] = imuraw.dvel[2];

    // Adjust data buffer index
    zuptindex_++;
    if (zuptindex_ >= imudatarate_) {
        zuptindex_ = 0;
    }
}

bool ZuptProcessor::IsZupt(const LcConfig& lcconfig, std::array<double, 6>* avg, double* quality) {
    std::array<double, 6> sum = {0};
    std::array<double, 6> std = {0};
    int gyrocnt = 0, acccnt = 0;
    static bool laststatus = false;

    avg->fill(0.0);
    for (int k = 0; k < imudatarate_; k++) {
        for (int i = 0; i < 6; i++) {
            (*avg)[i] += zuptbuff_[k * 6 + i];
        }
    }

    for (int i = 0; i < 6; i++) {
        (*avg)[i] *= lcconfig.imudatadt;
    }

    sum.fill(0.0);
    for (int k = 0; k < imudatarate_; k++) {
        for (int i = 0; i < 6; i++) {
            sum[i] += std::pow(zuptbuff_[k * 6 + i] - (*avg)[i], 2);
        }
    }

    for (int i = 0; i < 6; i++) {
        std[i] = std::sqrt(sum[i] * lcconfig.imudatadt);
    }

    // Only if both gyro and accelerometer have values below threshold, zero velocity condition is met
    if (std[0] < lcconfig.zuptthr[0]) gyrocnt++;
    if (std[1] < lcconfig.zuptthr[0]) gyrocnt++;
    if (std[2] < lcconfig.zuptthr[0]) gyrocnt++;
    if (std[3] < lcconfig.zuptthr[1]) acccnt++;
    if (std[4] < lcconfig.zuptthr[1]) acccnt++;
    if (std[5] < lcconfig.zuptthr[1]) acccnt++;

    // Adjust observation weight based on zero velocity condition, output variance adjustment
    if (gyrocnt > 0 && acccnt > 0) {
        *quality = 36.0 / ((gyrocnt + acccnt) * (gyrocnt + acccnt));
        if (laststatus) {
            return true;
        }
        laststatus = true;
    } else {
        laststatus = false;
    }
    return false;
}