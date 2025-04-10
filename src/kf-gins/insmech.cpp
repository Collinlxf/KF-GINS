/*
 * KF-GINS: An EKF-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Liqiang Wang
 *    Contact : wlq@whu.edu.cn
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

#include "common/earth.h"
#include "common/rotation.h"

#include "insmech.h"

void INSMech::insMech(const PVA &pvapre, PVA &pvacur, const IMU &imupre, const IMU &imucur) {

    // perform velocity update, position updata and attitude update in sequence, irreversible order
    // 依次进行速度更新、位置更新、姿态更新, 不可调换顺序
    velUpdate(pvapre, pvacur, imupre, imucur);
    posUpdate(pvapre, pvacur, imupre, imucur);
    attUpdate(pvapre, pvacur, imupre, imucur);
}

void INSMech::velUpdate(const PVA &pvapre, PVA &pvacur, const IMU &imupre, const IMU &imucur) {

    Eigen::Vector3d d_vfb, d_vfn, d_vgn, gl, midvel, midpos;
    Eigen::Vector3d temp1, temp2, temp3;
    Eigen::Matrix3d cnn, I33 = Eigen::Matrix3d::Identity();
    Eigen::Quaterniond qne, qee, qnn, qbb, q1, q2;

    // 计算地理参数，子午圈半径和卯酉圈半径，地球自转角速度投影到n系, n系相对于e系转动角速度投影到n系，重力值
    // calculate geographic parameters, Meridian and Mao unitary radii,
    // earth rotational angular velocity projected to n-frame,
    // rotational angular velocity of n-frame to e-frame projected to n-frame, and gravity
    Eigen::Vector2d rmrn = Earth::meridianPrimeVerticalRadius(pvapre.pos(0));
    Eigen::Vector3d wie_n, wen_n;
    // 地球自转轴与东向正交，所以地球自转投影到n系时，只有北向和地向的分量
    wie_n << WGS84_WIE * cos(pvapre.pos[0]), 0, -WGS84_WIE * sin(pvapre.pos[0]);
    // wen_n 表示导航坐标系相对于地球坐标系的旋转角速度
    // wen_n[0]：东向速度引起的北向角速度分量。
    // wen_n[1]：北向速度引起的东向角速度分量。
    // wen_n[2]：东向速度引起的地向角速度分量。
    /*
    1. 地向速度表示垂直于地表的运动，地向速度主要与高度变化有关，
        而高度变化不直接影响地表的水平运动（北向和东向速度）。因此，地向速度在惯性导航计算中通常是独立处理的。
    2. 地球的曲率和旋转效应会导致东向速度和北向速度相互影响。这是由于当载体在地表移动时，地球的曲率和自转会引起速度分量的相互耦合。
    */
    wen_n << pvapre.vel[1] / (rmrn[1] + pvapre.pos[2]), -pvapre.vel[0] / (rmrn[0] + pvapre.pos[2]),
        -pvapre.vel[1] * tan(pvapre.pos[0]) / (rmrn[1] + pvapre.pos[2]);
    double gravity = Earth::gravity(pvapre.pos);

    // 旋转效应和双子样划桨效应
    // rotational and sculling motion
    temp1 = imucur.dtheta.cross(imucur.dvel) / 2;
    temp2 = imupre.dtheta.cross(imucur.dvel) / 12;
    temp3 = imupre.dvel.cross(imucur.dtheta) / 12;

    // b系比力积分项
    // velocity increment due to the specific force
    d_vfb = imucur.dvel + temp1 + temp2 + temp3;

    // 比力积分项投影到n系
    // velocity increment dut to the specfic force projected to the n-frame
    temp1 = (wie_n + wen_n) * imucur.dt / 2;
    cnn   = I33 - Rotation::skewSymmetric(temp1);
    d_vfn = cnn * pvapre.att.cbn * d_vfb;

    // 计算重力/哥式积分项
    // velocity increment due to the gravity and Coriolis force
    gl << 0, 0, gravity;
    d_vgn = (gl - (2 * wie_n + wen_n).cross(pvapre.vel)) * imucur.dt;

    // 得到中间时刻速度
    // velocity at k-1/2
    midvel = pvapre.vel + (d_vfn + d_vgn) / 2;

    // 外推得到中间时刻位置
    // position extrapolation to k-1/2
    qnn = Rotation::rotvec2quaternion(temp1);
    temp2 << 0, 0, -WGS84_WIE * imucur.dt / 2;
    qee = Rotation::rotvec2quaternion(temp2);
    qne = Earth::qne(pvapre.pos);
    qne = (qee * qne * qnn).normalized();
    midpos[2] = pvapre.pos[2] - midvel[2] * imucur.dt / 2;
    midpos    = Earth::blh(qne, midpos[2]);

    // 重新计算中间时刻的rmrn, wie_e, wen_n
    // recompute rmrn, wie_n, and wen_n at k-1/2
    rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
    wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 重新计算n系下平均比力积分项
    // recompute d_vfn
    temp3 = (wie_n + wen_n) * imucur.dt / 2;
    cnn   = I33 - Rotation::skewSymmetric(temp3);
    d_vfn = cnn * pvapre.att.cbn * d_vfb;

    // 重新计算重力、哥式积分项
    // recompute d_vgn
    gl << 0, 0, Earth::gravity(midpos);
    d_vgn = (gl - (2 * wie_n + wen_n).cross(midvel)) * imucur.dt;

    // 速度更新完成
    // velocity update finish
    pvacur.vel = pvapre.vel + d_vfn + d_vgn;
}

void INSMech::posUpdate(const PVA &pvapre, PVA &pvacur, const IMU &imupre, const IMU &imucur) {

    Eigen::Vector3d temp1, temp2, midvel, midpos;
    Eigen::Quaterniond qne, qee, qnn;

    // 重新计算中间时刻的速度和位置
    // recompute velocity and position at k-1/2
    midvel = (pvacur.vel + pvapre.vel) / 2;
    midpos = pvapre.pos + Earth::DRi(pvapre.pos) * midvel * imucur.dt / 2;

    // 重新计算中间时刻地理参数
    // recompute rmrn, wie_n, wen_n at k-1/2
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
    wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 重新计算 k时刻到k-1时刻 n系旋转矢量
    // recompute n-frame rotation vector (n(k) with respect to n(k-1)-frame)
    temp1 = (wie_n + wen_n) * imucur.dt;
    /*
    1. 导航坐标系（n系）的旋转矢量表示导航坐标系相对于地球坐标系的旋转
    2. 导航坐标系的旋转是相对于前一个时刻的，因此表示从k时刻到k-1时刻的旋转
    */
    qnn   = Rotation::rotvec2quaternion(temp1);
    // e系转动等效旋转矢量 (k-1时刻k时刻，所以取负号)
    // e-frame rotation vector (e(k-1) with respect to e(k)-frame)
    /*
    1. 地球坐标系（e系）的旋转矢量表示地球的自转的影响。
    2. 地球自转是一个连续的过程，从k-1时刻到k时刻，地球自转了一个角度
    3. 因此，e系的旋转矢量表示从k-1时刻到k时刻地球自转的影响。
    */
    temp2 << 0, 0, -WGS84_WIE * imucur.dt;
    qee = Rotation::rotvec2quaternion(temp2);

    // 位置更新完成
    // position update finish
    // 前一个时刻n系到e系的旋转四元数
    qne           = Earth::qne(pvapre.pos);
    // 当前时刻n系到e系旋转四元数 = 两时刻e系旋转四元数 * 先前n系到e系旋转四元数 * 两时刻n系旋转四元数：
    qne           = (qee * qne * qnn).normalized();
    // 当前时刻高程 = 先前高程 - 高程方向速度 * 采样间隔（因为北东地，计算出的速度时地向的，所以减）：
    pvacur.pos[2] = pvapre.pos[2] - midvel[2] * imucur.dt;
    // 调用 blh() 根据 n 系到 e 系旋转四元数计算经纬度：
    pvacur.pos    = Earth::blh(qne, pvacur.pos[2]);
}

void INSMech::attUpdate(const PVA &pvapre, PVA &pvacur, const IMU &imupre, const IMU &imucur) {

    Eigen::Quaterniond qne_pre, qne_cur, qne_mid, qnn, qbb;
    Eigen::Vector3d temp1, midpos, midvel;

    // 重新计算中间时刻的速度和位置
    // recompute velocity and position at k-1/2
    midvel = (pvapre.vel + pvacur.vel) / 2;
    qne_pre   = Earth::qne(pvapre.pos);
    qne_cur   = Earth::qne(pvacur.pos);
    temp1     = Rotation::quaternion2vector((qne_cur.inverse() * qne_pre).normalized());
    qne_mid   = (qne_pre * Rotation::rotvec2quaternion(temp1 / 2).inverse()).normalized();
    midpos[2] = (pvacur.pos[2] + pvapre.pos[2]) / 2;
    midpos    = Earth::blh(qne_mid, midpos[2]);

    // 重新计算中间时刻地理参数
    // recompute rmrn, wie_n, wen_n at k-1/2
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
    wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 计算n系的旋转四元数 k-1时刻到k时刻变换
    // n-frame rotation vector (n(k-1) with respect to n(k)-frame)
    temp1 = -(wie_n + wen_n) * imucur.dt;
    qnn   = Rotation::rotvec2quaternion(temp1);

    // 计算b系旋转四元数 补偿二阶圆锥误差
    // b-frame rotation vector (b(k) with respect to b(k-1)-frame)
    // compensate the second-order coning correction term.
    temp1 = imucur.dtheta + imupre.dtheta.cross(imucur.dtheta) / 12;
    qbb   = Rotation::rotvec2quaternion(temp1);

    // 姿态更新完成
    // attitude update finish
    pvacur.att.qbn   = (qnn * pvapre.att.qbn * qbb).normalized();
    pvacur.att.cbn   = Rotation::quaternion2matrix(pvacur.att.qbn);
    pvacur.att.euler = Rotation::matrix2euler(pvacur.att.cbn);
    Eigen::Vector3d ibwb;/* b系下的陀螺仪观测量 */
    ibwb = imucur.dtheta/imucur.dt;
    Eigen::Vector3d temp22;
    temp22 = wie_n + wen_n;
    Eigen::Vector3d temp11;
    temp11 = pvacur.att.cbn *ibwb;
    pvacur.nbwn = temp11 - temp22;
}