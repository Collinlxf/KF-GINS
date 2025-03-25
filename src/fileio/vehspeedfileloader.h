#pragma once

#include "common/angle.h"
#include "common/types.h"
#include "fileloader.h"

class VehSpeedFileLoader : public FileLoader {

public:
    VehSpeedFileLoader() = delete;
    // columns最好跟文件的实际列数一致
    explicit VehSpeedFileLoader(const string &filename, int columns = 2) {
        open(filename, columns, FileLoader::TEXT);
    }

    const Veh_Speed &next() {
        data_ = load();

        // 第一列是时间
        veh_spd_.time = data_[0];
        veh_spd_.speed_veh = data_[1];

        // std::cout << std::fixed << std::setprecision(5) << __FILE__ << __LINE__
        // << "veh_spd_.time " << veh_spd_.time << std::endl;
        // std::cout << __FILE__ << __LINE__ << "veh_spd_.speed_veh : " << veh_spd_.speed_veh  << std::endl;
        return veh_spd_;
    }

private:
    Veh_Speed veh_spd_;
    vector<double> data_;
};
