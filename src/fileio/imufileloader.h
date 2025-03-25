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

#ifndef IMUFILELOADER_H
#define IMUFILELOADER_H

#include "common/types.h"
#include "fileloader.h"

class ImuFileLoader : public FileLoader {

public:
    ImuFileLoader() = delete;
    ImuFileLoader(const string &filename, int columns, int rate = 50) {
        open(filename, columns, FileLoader::TEXT);

        // 根据采样率计算时间间隔
        dt_ = 1.0 / (double) rate;

        // 初始化imu_.time的数据
        imu_.time = 0;
    }

    const IMU &next() {
        imu_pre_ = imu_;

        data_ = load();

        imu_.time = data_[0];
        // 第二列到第四列赋值给dtheta（角增量）
        memcpy(imu_.dtheta.data(), &data_[1], 3 * sizeof(double));
        // 第五列到第七列赋值给dvel（速度增量）
        memcpy(imu_.dvel.data(), &data_[4], 3 * sizeof(double));

        // 计算时间差 dt，若时间差小于0.1，则使用实际时间差；否则使用固定时间间隔 dt_
        double dt = imu_.time - imu_pre_.time;
        if (dt < 0.1) {
            imu_.dt = dt;
        } else {
            imu_.dt = dt_;
        }

        // 增量形式
        // 如果数据大于7列，第八列赋值给odovel（里程增量）
        if (columns_ > 7) {
            imu_.odovel = data_[7] * imu_.dt;
        }

        return imu_;
    }

    double starttime() {
        // 获取文件开始是时间

        double starttime;
        std::streampos sp = filefp_.tellg();

        filefp_.seekg(0, std::ios_base::beg);
        starttime = load().front();
        filefp_.seekg(sp, std::ios_base::beg);
        return starttime;
    }

    double endtime() {
        // 获取文件结束时间

        double endtime    = -1;
        std::streampos sp = filefp_.tellg();

        if (filetype_ == TEXT) {
            filefp_.seekg(-2, std::ios_base::end);
            char byte = 0;
            auto pos  = filefp_.tellg();
            do {
                pos -= 1;
                filefp_.seekg(pos);
                filefp_.read(&byte, 1);
            } while (byte != '\n');
        } else {
            filefp_.seekg(-columns_ * sizeof(double), std::ios_base::end);
        }
        endtime = load().front();
        filefp_.seekg(sp, std::ios_base::beg);
        return endtime;
    }

private:
    double dt_;

    IMU imu_, imu_pre_;
    vector<double> data_;
};

#endif // IMUFILELOADER_H
