/*
 * OB_GINS: An Optimization-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Hailiang Tang
 *    Contact : thl@whu.edu.cn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef GNSSFILELOADER_H
#define GNSSFILELOADER_H

#include "common/angle.h"
#include "common/types.h"
#include "fileloader.h"

class GnssFileLoader : public FileLoader {

public:
    GnssFileLoader() = delete;
    // columns最好跟文件的实际列数一致
    explicit GnssFileLoader(const string &filename, int columns = 10) {
        open(filename, columns, FileLoader::TEXT);
    }

    const GNSS &next() {
        data_ = load();

        // 第一列是时间
        gnss_.time = data_[0];
        // 第二列到第四列是BLH坐标的纬度、经度、高程
        memcpy(gnss_.blh.data(), &data_[1], 3 * sizeof(double));

        // 13列GNSS文件包含GNSS速度
        if (data_.size() == 7) {
            // 数据是7列，那么后面三列是标准差
            memcpy(gnss_.std.data(), &data_[4], 3 * sizeof(double));
        } else if (data_.size() == 10) {
            // 数据是8列，那么第5、6、7列是标准差，第八列是航向角
            memcpy(gnss_.std.data(), &data_[4], 3 * sizeof(double));
            gnss_.std[0] = 5.0;
            gnss_.std[1] = 5.0;
            gnss_.std[2] = 8.0;
            gnss_.yaw = data_[7];
            gnss_.speed_gps = data_[8];
            gnss_.sat_num = data_[9];
        } else {
            memcpy(gnss_.std.data(), &data_[7], 3 * sizeof(double));
        }
        gnss_.blh[0] *= D2R;
        gnss_.blh[1] *= D2R;
        // std::cout << __FILE__ << __LINE__ << "gnss_.std[0]: " << gnss_.std[0] << std::endl;
        // std::cout << __FILE__ << __LINE__ << "gnss_.std[1]: " << gnss_.std[1] << std::endl;
        // std::cout << __FILE__ << __LINE__ << "gnss_.std[2]: " << gnss_.std[2] << std::endl;
        return gnss_;
    }

private:
    GNSS gnss_;
    vector<double> data_;
};

#endif // GNSSFILELOADER_H
