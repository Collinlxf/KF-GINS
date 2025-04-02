#include "common/earth.h"
#include "common/rotation.h"

#include "gi_engine.h"
#include "insmech.h"

GIEngine::GIEngine(GINSOptions &options) {

    this->options_ = options;
    options_.print_options();
    timestamp_ = 0;

    // 初始化轮速计数据结构
    veh_speed_ = {0.0, 0.0};

    /*
    位置误差 (P_ID) - 3维
    速度误差 (V_ID) - 3维
    姿态误差 (PHI_ID) - 3维，使用phi角误差模型
    陀螺仪零偏误差 (BG_ID) - 3维
    加速度计零偏误差 (BA_ID) - 3维
    陀螺仪比例因子误差 (SG_ID) - 3维
    加速度计比例因子误差 (SA_ID) - 3维
    */

    // 设置协方差矩阵，系统噪声阵和系统误差状态矩阵大小
    // resize covariance matrix, system noise matrix, and system error state matrix
    // 协方差矩阵：描述了状态估计的不确定性。其对角线元素表示各状态变量的方差，非对角线元素表示状态变量之间的协方差。
    // 在预测步骤中，协方差矩阵。会根据系统模型进行传播，从而预测下一时刻的协方差。
    Cov_.resize(RANK, RANK);
    // 系统噪声矩阵，描述了系统过程中的随机扰动和不确定性。这个矩阵反映了系统模型的不完美性和外界干扰。
    Qc_.resize(NOISERANK, NOISERANK);
    // 系统误差状态矩阵：描述了状态估计与真实状态之间的误差。这个矩阵通常用于表示状态估计的微小偏差。
    // 会在EKF的预测和更新过程中不断调整，以改进状态估计的准确性。
    dx_.resize(RANK, 1);
    Cov_.setZero();
    Qc_.setZero();
    dx_.setZero();

    // 初始化系统噪声阵，系统噪声矩阵一旦设定，在后续ekf过程中不需要改变
    // initialize noise matrix
    // 系统噪声通常是标准差，需要平方处理
    auto imunoise                   = options_.imunoise;
    Qc_.block(ARW_ID, ARW_ID, 3, 3) = imunoise.gyr_arw.cwiseProduct(imunoise.gyr_arw).asDiagonal();
    // std::cout << "Qc_ Matrix:\n" << Qc_ << std::endl;
    Qc_.block(VRW_ID, VRW_ID, 3, 3) = imunoise.acc_vrw.cwiseProduct(imunoise.acc_vrw).asDiagonal();
    // std::cout << "Qc_ Matrix:\n" << Qc_ << std::endl;
    Qc_.block(BGSTD_ID, BGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrbias_std.cwiseProduct(imunoise.gyrbias_std).asDiagonal();
    // std::cout << "Qc_ Matrix:\n" << Qc_ << std::endl;
    Qc_.block(BASTD_ID, BASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accbias_std.cwiseProduct(imunoise.accbias_std).asDiagonal();
    // std::cout << "Qc_ Matrix:\n" << Qc_ << std::endl;
    Qc_.block(SGSTD_ID, SGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrscale_std.cwiseProduct(imunoise.gyrscale_std).asDiagonal();
    // std::cout << "Qc_ Matrix:\n" << Qc_ << std::endl;
    Qc_.block(SASTD_ID, SASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accscale_std.cwiseProduct(imunoise.accscale_std).asDiagonal();
    // std::cout << "Qc_ Matrix:\n" << Qc_ << std::endl;

    // 设置系统状态(位置、速度、姿态和IMU误差)初值和初始协方差
    // set initial state (position, velocity, attitude and IMU error) and covariance
    initialize(options_.initstate, options_.initstate_std);
}

void GIEngine::initialize(const NavState &initstate, const NavState &initstate_std) {

    // 初始化位置、速度、姿态
    // initialize position, velocity and attitude
    pvacur_.pos       = initstate.pos;
    pvacur_.vel       = initstate.vel;
    /* 
    汽车类的组合导航中pvacur_.att.euler中的欧拉角通常表示从导航坐标系到载体坐标系的旋转变化
    gps的yaw角可以直接用作欧拉角的yaw角，但是pitch和roll角需要根据载体坐标系的速度和加速度计算得到
    gps的yaw角可以直接用作欧拉角的yaw角的原因：
    1. 载体坐标系的x轴是指向前方的
    2. gps提供的yaw角是车辆前进方向（车头方向）与地理北向的夹角
    3. 欧拉角的yaw角是载体坐标系的x轴（车头方向）与导航坐标系的北向的夹角
    */
    pvacur_.att.euler = initstate.euler;
    pvacur_.att.cbn   = Rotation::euler2matrix(pvacur_.att.euler);
    pvacur_.att.qbn   = Rotation::euler2quaternion(pvacur_.att.euler);
    // 初始化IMU误差
    // initialize imu error
    imuerror_ = initstate.imuerror;

    // 给上一时刻状态赋同样的初值
    // set the same value to the previous state
    pvapre_ = pvacur_;

    // 初始化协方差
    // initialize covariance
    /*
    协方差矩阵中的对角线元素表示各个状态变量的方差，而非对角线元素表示状态变量之间的协方差。对角线元素和非对角线元素的含义如下：
    1. 对角线元素：表示单个状态变量的方差，即该变量的估计误差的平方。方差越大，表示该变量的不确定性越高。
    2. 非对角线元素：表示两个状态变量之间的协方差，即它们之间的相关性。如果两个变量之间没有相关性，协方差为零。
    初始化的时候为什么只初始化对角线元素：
    1. 初始简化：在系统初始化时，假设各个误差源之间是独立的，可以简化计算和分析。随着系统运行，
        滤波算法（如卡尔曼滤波）会根据实际数据更新协方差矩阵，包括非对角线元素。
    2. 无先验信息：在很多情况下，初始时没有先验信息表明各个状态变量之间存在显著的相关性。因此，默认各个变量之间没有关联。
    3. 递推更新：滤波器（如卡尔曼滤波器）在后续的状态估计过程中，会根据传感器数据和系统模型递推更新协方差矩阵，包括可能的相关性。
    */
    ImuError imuerror_std            = initstate_std.imuerror;
    Cov_.block(P_ID, P_ID, 3, 3)     = initstate_std.pos.cwiseProduct(initstate_std.pos).asDiagonal();
    Cov_.block(V_ID, V_ID, 3, 3)     = initstate_std.vel.cwiseProduct(initstate_std.vel).asDiagonal();
    Cov_.block(PHI_ID, PHI_ID, 3, 3) = initstate_std.euler.cwiseProduct(initstate_std.euler).asDiagonal();
    Cov_.block(BG_ID, BG_ID, 3, 3)   = imuerror_std.gyrbias.cwiseProduct(imuerror_std.gyrbias).asDiagonal();
    Cov_.block(BA_ID, BA_ID, 3, 3)   = imuerror_std.accbias.cwiseProduct(imuerror_std.accbias).asDiagonal();
    Cov_.block(SG_ID, SG_ID, 3, 3)   = imuerror_std.gyrscale.cwiseProduct(imuerror_std.gyrscale).asDiagonal();
    Cov_.block(SA_ID, SA_ID, 3, 3)   = imuerror_std.accscale.cwiseProduct(imuerror_std.accscale).asDiagonal();
}

void GIEngine::newImuProcess() {

    // 当前IMU时间作为系统当前状态时间,
    // set current IMU time as the current state time
    timestamp_ = imucur_.time;
    // 判断是否需要进行GNSS更新
    gnss_valid_ =  gnssdata_.sat_num > 5;
    std::cout << __FILE__ << __LINE__ << "sat_num: " << int(gnssdata_.sat_num) << std::endl;
    std::cout << __FILE__ << __LINE__ << "gnss_valid_: " << gnss_valid_ << std::endl;

    // 如果GNSS有效，则将更新时间设置为GNSS时间
    // set update time as the gnss time if gnssdata is valid
    double updatetime = gnssdata_.isvalid ? gnssdata_.time : -1;

    // 判断是否需要进行GNSS更新
    // determine if we should do GNSS update
    int res = isToUpdate(imupre_.time, imucur_.time, updatetime);
    // res = 0;

    if (res == 0) {
        // 只传播导航状态
        // only propagate navigation state
        std::cout << __FILE__ << __LINE__ << "res: " << res << std::endl;
        double gnss_yaw = Angle::deg2rad(gnssdata_.yaw);
        // pvacur_.att.euler[2] = gnss_yaw;
        // std::cout << __FILE__ << __LINE__ << "gnss_yaw: " << gnss_yaw << std::endl;
        insPropagation(imupre_, imucur_);
        // bool odd_valid = true;
        // if ((false == gnss_valid_ || (imucur_.time  > 641159.49100 && imucur_.time  < 641218.5)
        //     && (odd_valid == true && veh_speed_.time - imucur_.time < 0.03))) {
        //     std::cout << __FILE__ << __LINE__ << "ODONHCSpeedUpdate: "  << std::endl;
        //     // ODONHCSpeedUpdate();
        //     lcSpeedUpdate(pvacur_, veh_speed_.speed_veh);
        //     gnss_valid_ = false;
        //     odd_valid = false;
        // }

        bool is_veh_update = std::trunc(veh_speed_.time) != std::trunc(veh_speed_pre_.time);
        if (true == is_veh_update) {
            std::cout << __FILE__ << __LINE__ << "veh_speed_.time: " << veh_speed_.time << std::endl;
            std::cout << __FILE__ << __LINE__ << "veh_speed_pre_.time: " << veh_speed_pre_.time << std::endl;
        }
        if ((false == gnss_valid_ || (imucur_.time  > 641159.49100 && imucur_.time  < 641218.5)) && is_veh_update) {
            std::cout << __FILE__ << __LINE__ << "ODONHCSpeedUpdate: "  << std::endl;
            ODONHCSpeedUpdate();
            // lcSpeedUpdate(pvacur_, veh_speed_.speed_veh);
            gnss_valid_ = false;
        }
        if (false == gnss_valid_ || (imucur_.time  > 641159.49100 && imucur_.time  < 641218.5)) {
            gnss_valid_ = false;
        }

        // ODONHCSpeedUpdate();
        // std::cout << __FILE__ << __LINE__ << "pvacur_.rool: " << pvacur_.vel[0] << std::endl;
        // std::cout << __FILE__ << __LINE__ << "pvacur_.pitch: " << pvacur_.vel[1] << std::endl;
        // std::cout << __FILE__ << __LINE__ << "pvacur_.yaw: " << pvacur_.att.euler[2] << std::endl;
    } else if (res == 1) {
        // GNSS数据靠近上一历元，先对上一历元进行GNSS更新
        // gnssdata is near to the previous imudata, we should firstly do gnss update
        std::cout << __FILE__ << __LINE__ << "res: " << res << std::endl;
        // gnssUpdate进行 GNSS 量测更新
        // gnssUpdate(gnssdata_);
        // // stateFeedback进行系统状态反馈
        // stateFeedback();
        // ODONHCSpeedUpdate();

        pvapre_ = pvacur_;
        insPropagation(imupre_, imucur_);
    } else if (res == 2) {
        // GNSS数据靠近当前历元，先对当前IMU进行状态传播
        // gnssdata is near current imudata, we should firstly propagate navigation state
        std::cout << __FILE__ << __LINE__ << "res: " << res << std::endl;
        double gnss_yaw = Angle::deg2rad(gnssdata_.yaw);
        // std::cout << __FILE__ << __LINE__ << "gnss_yaw: " << gnss_yaw << std::endl;
        // pvacur_.att.euler[2] = gnss_yaw;
        // ODONHCSpeedUpdate();
        insPropagation(imupre_, imucur_);
        // std::cout << __FILE__ << __LINE__ << "pvacur_.rool: " << pvacur_.vel[0] << std::endl;
        // std::cout << __FILE__ << __LINE__ << "pvacur_.pitch: " << pvacur_.vel[1] << std::endl;
        // std::cout << __FILE__ << __LINE__ << "pvacur_.yaw: " << pvacur_.att.euler[2] << std::endl;
        // gnssUpdate(gnssdata_);
        // stateFeedback();
    } else {
        // GNSS数据在两个IMU数据之间(不靠近任何一个), 将当前IMU内插到整秒时刻
        // gnssdata is between the two imudata, we interpolate current imudata to gnss time
        std::cout << __FILE__ << __LINE__ << "res: " << res << std::endl;
        IMU midimu;
        imuInterpolate(imupre_, imucur_, updatetime, midimu);

        // 对前一半IMU进行状态传播
        // propagate navigation state for the first half imudata
        insPropagation(imupre_, midimu);

        // 整秒时刻进行GNSS更新，并反馈系统状态
        // do GNSS position update at the whole second and feedback system states
        // gnssUpdate(gnssdata_);
        // stateFeedback();
        // ODONHCSpeedUpdate();

        // 对后一半IMU进行状态传播
        // propagate navigation state for the second half imudata
        pvapre_ = pvacur_;
        insPropagation(midimu, imucur_);
    }

    if (true == gnss_valid_) {
        std::cout << __FILE__ << __LINE__ << "gnssUpdate: "  << std::endl;
        gnssUpdate(gnssdata_);
        // dx_.setZero();
        stateFeedback();
        dx_.setZero();
    }

    // 检查协方差矩阵对角线元素
    // check diagonal elements of current covariance matrix
    checkCov();

    // 更新上一时刻的状态和IMU数据
    // update system state and imudata at the previous epoch
    pvapre_ = pvacur_;
    imupre_ = imucur_;
}

void GIEngine::wheelSpeedUpdate() {
    // 利用轮速计数据进行速度更新
    double wheel_speed = veh_speed_.speed_veh;
    Eigen::Vector3d wheel_vel = {wheel_speed, 0, 0};

    // 计算轮速计速度测量新息
    Eigen::MatrixXd dz = wheel_vel - pvacur_.vel;

    // 构造轮速计速度观测矩阵
    Eigen::MatrixXd H_wheel;
    H_wheel.resize(3, Cov_.rows());
    H_wheel.setZero();
    H_wheel.block(0, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 速度观测噪声阵
    Eigen::MatrixXd R_wheel;
    // R_wheel = Eigen::Matrix3d::Identity() * options_.wheel_noise;
    R_wheel = Eigen::Matrix3d::Identity() * 0.1;

    // EKF更新协方差和误差状态
    EKFUpdate(dz, H_wheel, R_wheel);
}

void GIEngine::updateWheelSpeed(const Veh_Speed &veh_speed) {
    veh_speed_ = veh_speed;
}

void GIEngine::imuCompensate(IMU &imu) {

    // 补偿IMU零偏
    // compensate the imu bias
    imu.dtheta -= imuerror_.gyrbias * imu.dt;
    imu.dvel -= imuerror_.accbias * imu.dt;

    // 补偿IMU比例因子
    // compensate the imu scale
    Eigen::Vector3d gyrscale, accscale;
    gyrscale   = Eigen::Vector3d::Ones() + imuerror_.gyrscale;
    accscale   = Eigen::Vector3d::Ones() + imuerror_.accscale;
    // cwiseProduct：逐元素乘法。cwiseInverse：逐元素取倒数。
    imu.dtheta = imu.dtheta.cwiseProduct(gyrscale.cwiseInverse());
    imu.dvel   = imu.dvel.cwiseProduct(accscale.cwiseInverse());
}

void GIEngine::insPropagation(IMU &imupre, IMU &imucur) {

    bool zero_speed = false;
    // isZeroSpeed(imucur_, gnssdata_, &zero_speed);
    // if (true == zero_speed) {
    if (true == false) {
        // handleZeroSpeedCorrection(imucur_);
        //  std::cout << __FILE__ << __LINE__ <<"Cov_ Matrix Diagonal Elements:\n" << Cov_.diagonal() << std::endl;
        return;
    } else {
        // 对当前IMU数据(imucur)补偿误差, 上一IMU数据(imupre)已经补偿过了
        // compensate imu error to 'imucur', 'imupre' has been compensated
        imuCompensate(imucur);
        // IMU状态更新(机械编排算法)
        // update imustate(mechanization)
        INSMech::insMech(pvapre_, pvacur_, imupre, imucur);
        //  std::cout << __FILE__ << __LINE__ << "Cov_ Matrix Diagonal Elements:\n" << Cov_.diagonal() << std::endl;
    }

    // 系统噪声传播，姿态误差采用phi角误差模型
    // system noise propagate, phi-angle error model for attitude error
    Eigen::MatrixXd Phi, F, Qd, G;

    // 初始化Phi阵(状态转移矩阵)，F阵，Qd阵(传播噪声阵)，G阵(噪声驱动阵)
    // initialize Phi (state transition), F matrix, Qd(propagation noise) and G(noise driven) matrix
    /*
    矩阵	作用	            更新方式	                        更新时机
    Phi	状态转移预测	        基于IMU数据和地球模型动态计算	    每次IMU数据更新时
    F	线性化状态转移方程	    基于IMU数据和地球模型动态计算	    每次IMU数据更新时
    Qd	建模系统噪声统计特性	预设或基于IMU噪声参数动态调参	        初始化或调参阶段
    G	映射噪声到状态空间	    基于当前IMU姿态动态生成	            每次IMU数据更新时
    */

    /*
    Phi矩阵（状态转移矩阵）
    数学定义：离散状态方程的系数矩阵，描述误差状态如何随时间演化
    计算方法：通过一阶泰勒近似：Φ = I + F·dt
    物理意义：表示当前误差状态对下一时刻误差状态的影响
    作用：用于预测下一时刻的误差状态和协方差传播

    F矩阵（系统动力学矩阵）
    数学定义：连续时间误差状态方程的雅可比矩阵
    计算方法：基于INS误差传播方程进行导数计算
    物理意义：描述各误差状态之间的动态关系
    作用：用于构造离散状态转移矩阵Phi
    F矩阵中包含：
    位置误差相对于位置、速度的变化率
    速度误差相对于位置、速度、姿态、加速度计零偏和比例因子的变化率
    姿态误差相对于位置、速度、姿态、陀螺仪零偏和比例因子的变化率
    IMU误差参数（零偏、比例因子）的一阶高斯-马尔科夫过程模型

    Qd矩阵（离散过程噪声协方差矩阵）
    数学定义：系统噪声的协方差矩阵
    计算方法：Qd = G·Qc·G'·dt + 高阶修正项
    物理意义：表征系统不确定性的增长程度
    作用：用于协方差矩阵的预测更新

    G矩阵（噪声驱动矩阵）
    数学定义：将噪声向量映射到状态空间的矩阵
    计算方法：基于当前姿态和IMU误差模型构建
    物理意义：描述各种噪声源（ARW、VRW等）如何影响系统状态
    作用：用于计算系统噪声协方差矩阵Qd
    */

    Phi.resizeLike(Cov_);
    F.resizeLike(Cov_);
    Qd.resizeLike(Cov_);
    G.resize(RANK, NOISERANK);
    Phi.setIdentity();
    F.setZero();
    Qd.setZero();
    G.setZero();

    // 使用上一历元状态计算状态转移矩阵
    // compute state transition matrix using the previous state
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    double gravity;
    rmrn    = Earth::meridianPrimeVerticalRadius(pvapre_.pos[0]);
    gravity = Earth::gravity(pvapre_.pos);
    wie_n << WGS84_WIE * cos(pvapre_.pos[0]), 0, -WGS84_WIE * sin(pvapre_.pos[0]);
    wen_n << pvapre_.vel[1] / (rmrn[1] + pvapre_.pos[2]), -pvapre_.vel[0] / (rmrn[0] + pvapre_.pos[2]),
        -pvapre_.vel[1] * tan(pvapre_.pos[0]) / (rmrn[1] + pvapre_.pos[2]);

    Eigen::Matrix3d temp;
    Eigen::Vector3d accel, omega;
    double rmh, rnh;

    rmh   = rmrn[0] + pvapre_.pos[2];
    rnh   = rmrn[1] + pvapre_.pos[2];
    accel = imucur.dvel / imucur.dt;
    omega = imucur.dtheta / imucur.dt;

    // 位置误差
    // position error
    /*
    δṗ = F_pp · δp + F_pv · δv
    δp 是位置误差向量 [δL, δλ, δh]（纬度误差、经度误差、高度误差）
    δv 是速度误差向量 [δvn, δve, δvd]（北向、东向、下向速度误差）
    F_pp 是位置误差对位置误差变化率的影响矩阵
    F_pv 是速度误差对位置误差变化率的影响矩阵
    */
   /*
   数据示例
   Eigen::MatrixXd F = Eigen::MatrixXd::Zero(10, 10);
   Eigen::Matrix3d temp;
    temp << 10, 20, 30,
            40, 50, 60,
            70, 80, 90;
    // 将 temp 的值赋给 F 的特定位置
    F.block(P_ID, P_ID, 3, 3) = temp;

    // 将单位矩阵赋给 F 的另一个特定位置
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    输出：
    Matrix F:
    10 20 30  1  0  0  0  0  0  0
    40 50 60  0  1  0  0  0  0  0
    70 80 90  0  0  1  0  0  0  0
    0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0
   */
    temp.setZero();
    /*
    纬度误差变化率(δL̇)对位置的依赖：
    -vd/(Rm+h)向下速度越大，纬度误差变化越快，表示下沉越快会导致纬度计算偏差越大
    vn/(Rm+h)：沿子午圈方向（北向）运动时，高度误差会影响纬度计算
    */ 
    temp(0, 0)                = -pvapre_.vel[2] / rmh;
    temp(0, 2)                = pvapre_.vel[0] / rmh;
    /*
    经度误差变化率(δλ̇)对位置的依赖：
    ve*tan(L)/(Rn+h)：东向速度在高纬度地区（tan(L)大）导致的经度误差变化更大，这反映了经线汇聚效应
    -(vd+vn*tan(L))/(Rn+h)：既与下沉速度有关，也与北向速度在高纬度区域的影响有关
    ve/(Rn+h)：高度误差会影响东向运动引起的经度计算
    */
    temp(1, 0)                = pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1)                = -(pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2)                = pvapre_.vel[1] / rnh;
    F.block(P_ID, P_ID, 3, 3) = temp;
    // 速度误差对位置误差的直接影响：
    /*
    单位矩阵表示速度误差直接积分成位置误差，即：
    北向速度误差导致纬度误差增长
    东向速度误差导致经度误差增长
    垂直速度误差导致高度误差增长
    */
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 速度误差相对于位置、速度、姿态、加速度计零偏和比例因子的变化率
    // velocity error
    temp.setZero();
    temp(0, 0) = -2 * pvapre_.vel[1] * WGS84_WIE * cos(pvapre_.pos[0]) / rmh -
                 pow(pvapre_.vel[1], 2) / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(0, 2) = pvapre_.vel[0] * pvapre_.vel[2] / rmh / rmh - pow(pvapre_.vel[1], 2) * tan(pvapre_.pos[0]) / rnh / rnh;
    temp(1, 0) = 2 * WGS84_WIE * (pvapre_.vel[0] * cos(pvapre_.pos[0]) - pvapre_.vel[2] * sin(pvapre_.pos[0])) / rmh +
                 pvapre_.vel[0] * pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(1, 2) = (pvapre_.vel[1] * pvapre_.vel[2] + pvapre_.vel[0] * pvapre_.vel[1] * tan(pvapre_.pos[0])) / rnh / rnh;
    temp(2, 0) = 2 * WGS84_WIE * pvapre_.vel[1] * sin(pvapre_.pos[0]) / rmh;
    temp(2, 2) = -pow(pvapre_.vel[1], 2) / rnh / rnh - pow(pvapre_.vel[0], 2) / rmh / rmh +
                 2 * gravity / (sqrt(rmrn[0] * rmrn[1]) + pvapre_.pos[2]);
    F.block(V_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 0)                  = pvapre_.vel[2] / rmh;
    temp(0, 1)                  = -2 * (WGS84_WIE * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh);
    temp(0, 2)                  = pvapre_.vel[0] / rmh;
    temp(1, 0)                  = 2 * WGS84_WIE * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1)                  = (pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2)                  = 2 * WGS84_WIE * cos(pvapre_.pos[0]) + pvapre_.vel[1] / rnh;
    temp(2, 0)                  = -2 * pvapre_.vel[0] / rmh;
    temp(2, 1)                  = -2 * (WGS84_WIE * cos(pvapre_.pos(0)) + pvapre_.vel[1] / rnh);
    F.block(V_ID, V_ID, 3, 3)   = temp;
    F.block(V_ID, PHI_ID, 3, 3) = Rotation::skewSymmetric(pvapre_.att.cbn * accel);
    F.block(V_ID, BA_ID, 3, 3)  = pvapre_.att.cbn;
    F.block(V_ID, SA_ID, 3, 3)  = pvapre_.att.cbn * (accel.asDiagonal());

    // 姿态误差相对于位置、速度、姿态、陀螺仪零偏和比例因子的变化率
    // attitude error
    temp.setZero();
    temp(0, 0) = -WGS84_WIE * sin(pvapre_.pos[0]) / rmh;
    temp(0, 2) = pvapre_.vel[1] / rnh / rnh;
    temp(1, 2) = -pvapre_.vel[0] / rmh / rmh;
    temp(2, 0) = -WGS84_WIE * cos(pvapre_.pos[0]) / rmh - pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(2, 2) = -pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh / rnh;
    F.block(PHI_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 1)                    = 1 / rnh;
    temp(1, 0)                    = -1 / rmh;
    temp(2, 1)                    = -tan(pvapre_.pos[0]) / rnh;
    F.block(PHI_ID, V_ID, 3, 3)   = temp;
    F.block(PHI_ID, PHI_ID, 3, 3) = -Rotation::skewSymmetric(wie_n + wen_n);
    F.block(PHI_ID, BG_ID, 3, 3)  = -pvapre_.att.cbn;
    F.block(PHI_ID, SG_ID, 3, 3)  = -pvapre_.att.cbn * (omega.asDiagonal());

    // IMU零偏误差和比例因子误差，建模成一阶高斯-马尔科夫过程
    // imu bias error and scale error, modeled as the first-order Gauss-Markov process
    F.block(BG_ID, BG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(BA_ID, BA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SG_ID, SG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SA_ID, SA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();

    // 系统噪声驱动矩阵
    // system noise driven matrix
    // G 只表示瞬时的噪声影响，不需要提前乘以 Δt
    // 噪声驱动阵里面只需要旋转矩阵举可以，最终产生的误差是噪声驱动阵（G）乘以系统噪声协方差矩阵（Qc_）
    // 位置误差 (P_ID) 没有直接受噪声影响，而是由速度误差积分得到，因此 G 矩阵没有针对 P_ID 进行单独赋值
    /*
    在扩展卡尔曼滤波（EKF）或无迹卡尔曼滤波（UKF）中，G 矩阵用于表示噪声 𝑤 对系统状态 𝑥的影响：
    𝑑𝑥 = 𝐹𝑥 + 𝐺𝑤
    其中：
        𝐺是 噪声驱动矩阵，描述过程噪声 𝑤如何影响状态变量 𝑥
        𝑤可能包括陀螺仪噪声、加速度计噪声、偏置漂移噪声等。
    */
    G.block(V_ID, VRW_ID, 3, 3)    = pvapre_.att.cbn;
    G.block(PHI_ID, ARW_ID, 3, 3)  = pvapre_.att.cbn;
    G.block(BG_ID, BGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(BA_ID, BASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SG_ID, SGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SA_ID, SASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 状态转移矩阵
    // compute the state transition matrix
    /*
    矩阵指数的泰勒展开，并截取到一阶项来近似计算状态转移矩阵。
    适用情况：
    1.Δt 很小时（如 IMU 数据的高频采样情况）。（如果Δt 比较大，那么需要用更高阶的泰勒展开）
    2.F 变化不剧烈的情况下
    */ 
    Phi.setIdentity();
    Phi = Phi + F * imucur.dt;
    // TODO 二阶泰勒展开
    // Phi = Phi + F * imucur.dt + 0.5 * F * F * imucur.dt * imucur.dt;

    // 计算系统传播噪声
    // compute system propagation noise
    // 𝑄𝑑作用于协方差 𝑃控制不确定性的传播。
    //  零阶保持近似
    Qd = G * Qc_ * G.transpose() * imucur.dt;
    // 修正过程噪声离散化误差，提高数值精度。van Loan 方法的近似公式，用于更精确地计算状态转移对噪声的影响。
    Qd = (Phi * Qd * Phi.transpose() + Qd) / 2;

    // EKF预测传播系统协方差和系统误差状态
    // do EKF predict to propagate covariance and error state
    // std::cout<<__FILE__<<","<<__LINE__<<", "<<" before EKFPredict: "<<std::endl;
    // std::cout<<__FILE__<<","<<__LINE__<<", "<<" before EKFPredict: dx_: \n"<<dx_<<std::endl;
    // std::cout<<__FILE__<<","<<__LINE__<<", "<<" before EKFPredict: cov_: \n"<<Cov_<<std::endl;
    EKFPredict(Phi, Qd);
    // std::cout<<__FILE__<<","<<__LINE__<<", "<<" before EKFPredict: "<<std::endl;
    // std::cout<<__FILE__<<","<<__LINE__<<", "<<" before EKFPredict: dx_: \n"<<dx_<<std::endl;
    // std::cout<<__FILE__<<","<<__LINE__<<", "<<" before EKFPredict: cov_: \n"<<Cov_<<std::endl;
}

void GIEngine::gnssUpdate(GNSS &gnssdata) {

    // IMU位置转到GNSS天线相位中心位置
    // convert IMU position to GNSS antenna phase center position
    Eigen::Vector3d antenna_pos;
    Eigen::Matrix3d Dr, Dr_inv;
    Dr_inv      = Earth::DRi(pvacur_.pos);
    Dr          = Earth::DR(pvacur_.pos);
    antenna_pos = pvacur_.pos + Dr_inv * pvacur_.att.cbn * options_.antlever;

    // GNSS位置测量新息
    // compute GNSS position innovation
    Eigen::MatrixXd dz;
    dz = Dr * (antenna_pos - gnssdata.blh);

    // 构造GNSS位置观测矩阵
    // construct GNSS position measurement matrix
    Eigen::MatrixXd H_gnsspos;
    /*
    观测矩阵 H 的列数需与协方差矩阵 Cov_ 的行数一致。这是因为矩阵乘法 H * Cov_ 要求 H 的列数等于 Cov_ 的行数。
    */
    H_gnsspos.resize(3, Cov_.rows());
    H_gnsspos.setZero();
    H_gnsspos.block(0, P_ID, 3, 3)   = Eigen::Matrix3d::Identity();
    H_gnsspos.block(0, PHI_ID, 3, 3) = Rotation::skewSymmetric(pvacur_.att.cbn * options_.antlever);

    // 位置观测噪声阵
    // construct measurement noise matrix
    Eigen::MatrixXd R_gnsspos;
    R_gnsspos = gnssdata.std.cwiseProduct(gnssdata.std).asDiagonal();

    // EKF更新协方差和误差状态
    // do EKF update to update covariance and error state
    EKFUpdate(dz, H_gnsspos, R_gnsspos);

    // GNSS更新之后设置为不可用
    // Set GNSS invalid after update
    gnssdata.isvalid = false;
}

int GIEngine::isToUpdate(double imutime1, double imutime2, double updatetime) const {

    if (abs(imutime1 - updatetime) < TIME_ALIGN_ERR) {
        // 更新时间靠近imutime1
        // updatetime is near to imutime1
        return 1;
    } else if (abs(imutime2 - updatetime) <= TIME_ALIGN_ERR) {
        // 更新时间靠近imutime2
        // updatetime is near to imutime2
        return 2;
    } else if (imutime1 < updatetime && updatetime < imutime2) {
        // 更新时间在imutime1和imutime2之间, 但不靠近任何一个
        // updatetime is between imutime1 and imutime2, but not near to either
        return 3;
    } else {
        // 更新时间不在imutimt1和imutime2之间，且不靠近任何一个
        // updatetime is not bewteen imutime1 and imutime2, and not near to either.
        return 0;
    }
}

void GIEngine::EKFPredict(Eigen::MatrixXd &Phi, Eigen::MatrixXd &Qd) {

    assert(Phi.rows() == Cov_.rows());
    assert(Qd.rows() == Cov_.rows());

    // 传播系统协方差和误差状态
    // propagate system covariance and error state
    // Pk+1 =Φk * Pk * ΦkT + Qd
    Cov_ = Phi * Cov_ * Phi.transpose() + Qd;
    // std::cout << __FILE__ << __LINE__ <<"Cov_ Matrix Diagonal Elements:\n" << Cov_.diagonal() << std::endl;
    dx_  = Phi * dx_;
}

void GIEngine::EKFUpdate(Eigen::MatrixXd &dz, Eigen::MatrixXd &H, Eigen::MatrixXd &R) {

    assert(H.cols() == Cov_.rows());
    assert(dz.rows() == H.rows());
    assert(dz.rows() == R.rows());
    assert(dz.cols() == 1);

    // 计算Kalman增益
    // Compute Kalman Gain
    auto temp         = H * Cov_ * H.transpose() + R;
    Eigen::MatrixXd K = Cov_ * H.transpose() * temp.inverse();

    // 更新系统误差状态和协方差
    // update system error state and covariance
    Eigen::MatrixXd I;
    I.resizeLike(Cov_);
    I.setIdentity();
    I = I - K * H;
    // 如果每次更新后都进行状态反馈，则更新前dx_一直为0，下式可以简化为：dx_ = K * dz;
    // if state feedback is performed after every update, dx_ is always zero before the update
    // the following formula can be simplified as : dx_ = K * dz;
    dx_  = dx_ + K * (dz - H * dx_);
    Cov_ = I * Cov_ * I.transpose() + K * R * K.transpose();
    // std::cout << __FILE__ << __LINE__ <<"Cov_ Matrix Diagonal Elements:\n" << Cov_.diagonal() << std::endl;
}

void GIEngine::stateFeedback() {

    Eigen::Vector3d vectemp;

    // 位置误差反馈
    // posisiton error feedback
    /*
    P_ID：起始的行索引
    0,起始的列索引
    3,要复制的行数
    1,要复制的列数
    */
    /*
    位置和速度误差修正：（通过减法实现）
    IMU零偏和比例因子误差修正：（通过加法实现）
    */
    Eigen::Vector3d delta_r = dx_.block(P_ID, 0, 3, 1);
    Eigen::Matrix3d Dr_inv  = Earth::DRi(pvacur_.pos);
    pvacur_.pos -= Dr_inv * delta_r;

    // 速度误差反馈
    // velocity error feedback
    vectemp = dx_.block(V_ID, 0, 3, 1);
    pvacur_.vel -= vectemp;
    if (false == gnss_valid_ || (imucur_.time  > 641159.49100 && imucur_.time  < 641218.5)) {
        // pvacur_.vel(1) = 0.0;
    }


    // 姿态误差反馈
    // attitude error feedback
    vectemp                = dx_.block(PHI_ID, 0, 3, 1);
    Eigen::Quaterniond qpn = Rotation::rotvec2quaternion(vectemp);
    pvacur_.att.qbn        = qpn * pvacur_.att.qbn;
    pvacur_.att.cbn        = Rotation::quaternion2matrix(pvacur_.att.qbn);
    pvacur_.att.euler      = Rotation::matrix2euler(pvacur_.att.cbn);

    // IMU零偏误差反馈
    // IMU bias error feedback
    vectemp = dx_.block(BG_ID, 0, 3, 1);
    imuerror_.gyrbias += vectemp;
    vectemp = dx_.block(BA_ID, 0, 3, 1);
    imuerror_.accbias += vectemp;

    // IMU比例因子误差反馈
    // IMU sacle error feedback
    vectemp = dx_.block(SG_ID, 0, 3, 1);
    imuerror_.gyrscale += vectemp;
    vectemp = dx_.block(SA_ID, 0, 3, 1);
    imuerror_.accscale += vectemp;

    /* 里程计比例因子 */
    // imuerror_.odoscale += dx_.block(OS_ID, 0, 1, 1)(0);

    // 误差状态反馈到系统状态后,将误差状态清零
    // set 'dx' to zero after feedback error state to system state
    // dx_.setZero();
}

void GIEngine::ODONHCUpdateFeedback(Eigen::MatrixXd &dz, Eigen::MatrixXd &H, Eigen::MatrixXd &R) {
    // dx_.setZero();
    EKFUpdate(dz, H, R);
    stateFeedback();
    dx_.setZero();
    return;
}

NavState GIEngine::getNavState() {

    NavState state;

    state.pos      = pvacur_.pos;
    state.vel      = pvacur_.vel;
    state.euler    = pvacur_.att.euler;
    state.imuerror = imuerror_;

    return state;
}

// 零速检测函数
bool GIEngine::isZeroSpeed(const IMU &imu, const GNSS &gnss, bool *iszero) {
    // 加速度幅值检测
    double acc_mag = imu.dvel.norm();  // 注意：此处使用dvel（已补偿速度）
    // std::cout << "acc_mag: " << acc_mag << std::endl;
    // if (acc_mag < 0.15) {
    //     std::cout << "acc_mag: " << acc_mag << std::endl;
    //     *iszero = true;
    //     return iszero;
    // }

    // // 陀螺仪角速度检测
    // double gyro_mag = imu.dtheta.norm();
    // std::cout << "gyro_mag: " << gyro_mag << std::endl;
    // if (gyro_mag < 0.01)  {
    //     std::cout << "gyro_mag: " << gyro_mag << std::endl;
    //     *iszero = true;
    //     return iszero;
    // }
    // 轮速计检测（可选）

    // std::cout << "gnss.speed_gps: " << gnss.speed_gps << std::endl;
    if (gnss.speed_gps < 0.2) {
        std::cout << "gnss.speed_gps: " << gnss.speed_gps << std::endl;
        *iszero = true;
        return iszero;
    }
    // GPS信号无效
    return false;
}

// 零速修正处理
void GIEngine::handleZeroSpeedCorrection(const IMU &imu) {
    // 冻结位置和速度误差：将协方差矩阵对角线元素置零
    double epsilon = 1e-6;  // 很小的正数，用于确保协方差矩阵的正定性
     Cov_.block(P_ID, P_ID, 3, 3) = Eigen::Matrix3d::Identity() * epsilon;
    Cov_.block(V_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity() * epsilon;

    // 强制补偿IMU零偏和比例因子
    imuerror_.gyrbias.setZero();    // 清零陀螺仪零偏
    imuerror_.accbias.setZero();    // 清零加速度计零偏
    imuerror_.gyrscale.setZero();   // 清零陀螺仪比例因子
    imuerror_.accscale.setZero();   // 清零加速度计比例因子

    // 更新当前状态
    pvacur_.vel.setZero();  // 速度置零
    pvacur_.pos = pvapre_.pos;  // 位置保持前一状态
}

bool GIEngine::ODONHCSpeedUpdate() {
    // 1. 里程计更新准备：计算 odovel
    // if (true == options_.lcconfig.id_use_odo || true == options_.lcconfig.is_use_zupt) {
    //     // if (lcdata->imudata.odosechalfnum == 0) {
    //     if (0 == 0) {
    //         imucur_.odovel = lcdata->imudata.odofirhalfvel;
    //     // } else if (lcdata->imudata.odofirhalfnum == 0) {
    //     } else if (1 == 0) {
    //         lcdata->imudata.odovel = lcdata->imudata.odosechalfvel;
    //     } else {
    //         lcdata->imudata.odovel = 1.5 * lcdata->imudata.odosechalfvel / lcdata->imudata.odosechalfnum -
    //                                  0.5 * lcdata->imudata.odofirhalfvel / lcdata->imudata.odofirhalfnum;
    //     }
    //     // 清零累积变量
    //     lcdata->imudata.odosechalfvel = 0;
    //     lcdata->imudata.odosechalfnum = 0;
    //     lcdata->imudata.odofirhalfvel = 0;
    //     lcdata->imudata.odofirhalfnum = 0;
    // }
    imucur_.odovel[0] = veh_speed_.speed_veh;
    std::cout << __FILE__ << __LINE__ << ", " << "imucur_.odovel: " << imucur_.odovel << std::endl;
    // 2. 零速修正判断
    double zuptquality = 1.0;
    if (true == options_.lcconfig.is_use_zupt) {
        if (true == options_.lcconfig.is_use_zupt) {
            // 里程计零速判断：若 |odovel| < 0.0001，则认为处于零速状态
            lcdata_.is_zupt = (std::abs(imucur_.odovel[0]) < 0.0001);
        } else {
            // 这里可以调用更复杂的零速检测算法，简化处理为 false
            lcdata_.is_zupt = false;
        }
    }

    // 3. 构造观测矩阵、测量值及噪声矩阵（使用 Eigen）
    Eigen::MatrixXd dz, H, R;
    // int stateDim = lcdata->rank;  // 状态向量维数
    int stateDim = 21;  // 状态向量维数

    // if (true == lcdata_.is_zupt) {
    if (true == false) {
        // 如果处于零速状态，则构造 3 维速度观测
        dz.resize(3, 1);
        // 直接使用当前测量到的速度
        dz = pvacur_.vel;  // 3×1向量

        // 构造观测矩阵 H：假设状态中速度误差位于 V_ID 开始的 3 个元素
        H = Eigen::MatrixXd::Zero(3, stateDim);
        H.block(0, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

        // 构造观测噪声矩阵 R：3×3 对角矩阵，噪声方差 = obsstd[2]^2 * zuptquality
        // double noise = lcconfig->obsstd(2) * lcconfig->obsstd(2) * zuptquality;
        double noise = 0.1 * 0.1 * zuptquality;
        R = noise * Eigen::MatrixXd::Identity(3, 3);
    // } else if (lcconfig->isusenhc) {
    } else if (1) {
        /*
        // NHC 速度观测：在此简化为只观测 x 轴速度误差
        dz.resize(2, 1);
        Eigen::Vector3d wie_n;
        wie_n << WGS84_WIE * cos(pvapre_.pos[0]), 0, -WGS84_WIE * sin(pvapre_.pos[0]);
        // 将自车坐标系前向速度转换到导航坐标系
        Eigen::Vector3d vel_car_to_nav = pvapre_.att.cbn * imucur_.odovel;
        // dz(0, 0) = pvacur_.vel(0) -  (vel_car_to_nav(0) + wie_n(0) * imucur_.dt / 2);
        dz(0, 0) = pvapre_.vel(0) -  (vel_car_to_nav(0) + wie_n(0) * 0.02 / 2);
        std::cout << __FILE__ << __LINE__ << ", " << "pvacur_.vel(0): " << pvacur_.vel(0) << std::endl;
        std::cout << __FILE__ << __LINE__ << ", " << "vel_car_to_nav(0): " << vel_car_to_nav(0) << std::endl;
        std::cout << __FILE__ << __LINE__ << ", " << "wie_n(0): " << wie_n(0) << std::endl;
        std::cout << __FILE__ << __LINE__ << ", " << "dz: " << dz.transpose() << std::endl;

        Eigen::Matrix3d C_nb = pvapre_.att.cbn.transpose();  // 从导航到车体的旋转矩阵
        Eigen::Vector3d vel_body = C_nb * pvapre_.vel;
        // 侧向速度残差（期望值为0）
        double lateral_residual = 0.0 - vel_body(1);
        dz(1, 0) = lateral_residual;

        H = Eigen::MatrixXd::Zero(2, stateDim);
        // 观测矩阵 H 在速度状态中 x 方向的系数置 1（假设 x 轴速度对应状态向量中 V_ID 下标处）
        H(0, V_ID) = 1.0;
        double psi = -pvapre_.att.euler(2);
        // 根据车辆坐标系和导航系的转换关系：
        // v_y_body ≈ -sin(psi)*v_x_nav + cos(psi)*v_y_nav
        H(1, V_ID)     = -sin(psi);
        H(1, V_ID + 1) =  cos(psi);

        // double noise = lcconfig->obsstd(0) * lcconfig->obsstd(0);
        // double noise = 0.01 * 0.01;
        // R = noise * Eigen::MatrixXd::Identity(1, 1);
        // 3. 构造测量噪声矩阵 R (2x2)
        // 正向速度噪声与侧向约束噪声可以根据实验标定进行调整
        Eigen::MatrixXd R = Eigen::MatrixXd::Zero(2, 2);
        R(0, 0) = 0.1 * 0.1;   // 正向速度噪声
        R(1, 1) = 0.1 * 0.1;   // 侧向约束噪声
        */

        /*
        
        // 1. 构造 NHC 速度观测（无杆臂，b系与v系相同）
        // 此时，转换矩阵直接由导航到车体的旋转矩阵给出，即：
        Eigen::Matrix3d C_nb = pvapre_.att.cbn.transpose();  // 从导航到车体的旋转矩阵

        // 2. 计算速度转换
        // 将导航系下的速度转换到车体坐标系（这里得到的就是车辆在车体下的速度，包含正向和侧向分量）
        Eigen::Vector3d vel_body = C_nb * pvapre_.vel;

        // 此外，正向速度观测中还加入了地球旋转速率的影响（wie_n），保持不变
        Eigen::Vector3d wie_n;
        wie_n << WGS84_WIE * cos(pvapre_.pos[0]), 0, -WGS84_WIE * sin(pvapre_.pos[0]);

        // 3. 构造测量残差 dz
        // dz 的第一行为正向（x 轴）速度误差，第二行为侧向（车体 y 轴）速度误差（非完整约束要求侧向速度接近零）
        dz.resize(18, 1);
        dz.setZero();
        {
            // 将车体前向速度转换到导航系（注意：如果 imucur_.odovel 已经在车体系，则先转换到导航系）
            // 此处假设 imucur_.odovel 为车体系测量值
            Eigen::Vector3d vel_car_to_nav = pvapre_.att.cbn * imucur_.odovel;
            // 正向速度残差：期望值为 vel_car_to_nav(0) 加上地球旋转影响（近似补偿一半 dt 影响）
            // 注意 dt 这里固定使用 0.02（整秒或采样间隔，根据实际情况调整）
            dz(3, 0) = pvapre_.vel(0) - (vel_car_to_nav(0) + wie_n(0) * 0.02 / 2);

            // 侧向速度残差：期望值为 0（非完整约束，车体侧向速度应为 0）
            dz(4, 0) = 0.0 - vel_body(1);

            std::cout << __FILE__ << __LINE__ << ", " << "pvacur_.vel(0): " << pvacur_.vel(0) << std::endl;
            std::cout << __FILE__ << __LINE__ << ", " << "vel_car_to_nav(0): " << vel_car_to_nav(0) << std::endl;
            std::cout << __FILE__ << __LINE__ << ", " << "wie_n(0): " << wie_n(0) << std::endl;
            std::cout << __FILE__ << __LINE__ << ", " << "dz: " << dz.transpose() << std::endl;
        }

        // 4. 构造观测矩阵 H
        // 状态向量中假设速度分量在下标 V_ID（正向）、V_ID+1（横向）、V_ID+2（垂向）
        // 这里构造两行观测方程：
        // 第1行：正向速度误差，对应于 v_north 分量（或者直接使用正向速度分量，本例中直接置 1）
        // 第2行：侧向非完整约束，根据车辆坐标系与导航系转换关系确定
        H = Eigen::MatrixXd::Zero(18, stateDim);

        // 第一行：仅对正向速度（假设 x 轴）直接观测
        H(0, V_ID) = 1.0;

        // 第二行：由车辆坐标系转换关系得，假设航向角 psi 定义为：正北顺时针为正
        // 此处计算航向角并取反（根据实际定义可能需要调整符号）
        double psi = -pvapre_.att.euler(2);
        // 车辆侧向速度与导航系下速度的线性化关系：
        // v_y_body ≈ -sin(psi)*v_north + cos(psi)*v_east
        H(1, V_ID)     = -sin(psi);
        H(1, V_ID + 1) =  cos(psi);

        // 5. 构造测量噪声矩阵 R (2x2)
        // 根据实验标定设置噪声水平（此处均设置为 0.1m/s 的标准差）
        Eigen::MatrixXd R = Eigen::MatrixXd::Zero(18, 18);
        R(V_ID, V_ID) = 0.1 * 0.1;   // 正向速度噪声
        R(V_ID + 1, V_ID + 1) = 0.1 * 0.1;   // 侧向约束噪声

        */

        // 1. 计算从导航到车体的旋转矩阵
        Eigen::Matrix3d C_nb = pvapre_.att.cbn.transpose();  // C_nb = (C_bn)^T

        // 2. 将导航系下的速度转换到车体坐标系
        Eigen::Vector3d b_pva = C_nb * pvapre_.vel;  // 预测的车体速度

        // 3. 计算姿态敏感性项：apsi = C_nb * skew(velocity)
        // 这里的 skew(velocity) 为速度的反对称矩阵，描述小姿态误差对测量的影响
        Eigen::Matrix3d vel_skew = Rotation::skewSymmetric(pvapre_.vel);
        Eigen::Matrix3d apsi = C_nb * vel_skew; 

        // 4. 构造测量残差 dz (3×1)
        // 对于正向分量，用轮速计测得的 odo_speed 补偿；侧向与垂向分量期望为零
        Eigen::MatrixXd dz(3, 1);
        dz(0, 0) = b_pva(0) - imucur_.odovel(0);  // 前向速度残差
        dz(1, 0) = b_pva(1);              // 侧向残差
        dz(2, 0) = b_pva(2);              // 垂向残差

        // 5. 构造观测矩阵 H (3×stateDim)
        // 初始化 H 为零矩阵
        H = Eigen::MatrixXd::Zero(3, stateDim);

        // 5.1 对于速度状态 (状态下标 V_ID 到 V_ID+2)
        // 若需要调整测量顺序（例如将预测的 forward 对应为第一行），可以对 C_nb 进行重排
        // 例如，假设当前 C_nb 的行顺序为 [row0, row1, row2] 对应车体的 [右, 下, 前]，
        // 但我们希望按照 [前, 右, 下] 的顺序进行测量，则需要重排：
        Eigen::Matrix3d C_nb_fill;
        C_nb_fill.row(0) = C_nb.row(2);  // 前向分量
        C_nb_fill.row(1) = C_nb.row(0);  // 右向分量
        C_nb_fill.row(2) = C_nb.row(1);  // 下向分量
        // 将此部分赋值到 H 对应位置
        H.block(0, V_ID, 3, 3) = C_nb_fill;

        // 5.2 对于姿态状态 (状态下标 PHI_ID 到 PHI_ID+2)
        // 计算姿态部分的雅可比矩阵：-C_nb * skew(v)
        // 同样，对 apsi 进行相应的重排：
        Eigen::Matrix3d apsi_fill;
        apsi_fill.row(0) = -apsi.row(2);
        apsi_fill.row(1) = -apsi.row(0);
        apsi_fill.row(2) = -apsi.row(1);
        H.block(0, PHI_ID, 3, 3) = apsi_fill;

        // 6. 构造测量噪声矩阵 R (3×3)
        // 例如，前向测量噪声与侧向、垂向噪声分别设定为 0.1 m/s 的标准差
        Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);
        R(0, 0) = 0.1 * 0.1;
        R(1, 1) = 0.1 * 0.1;
        R(2, 2) = 0.1 * 0.1;

        // // 7. 调用 EKF 更新及状态反馈
        // EKFUpdate(dz, H, R);
        // stateFeedback();
        // dx_.setZero();

    } else {
        // 如果既不零速也不使用 NHC，则不做更新
        return true;
    }

    EKFUpdate(dz, H, R);
    stateFeedback();
    dx_.setZero();

    // // 4. 调用 EKF 更新函数进行协方差和状态更新
    // ODONHCUpdateFeedback(dz, H, R);

    return true;
}

bool GIEngine::lcSpeedUpdate(PVA pvacur, double odo_speed){
    /* 计算nbwb */  
    Eigen::Matrix3d cnb = pvacur.att.cbn.transpose();
    Eigen::Vector3d tmp3 = cnb*pvacur.nbwn;
  
    Eigen::Vector3d odolever(0,0,0);
    /* 里程计杆臂带来的速度差异 */
    tmp3 = {0,0,0};
  
    /* 在速度改变之前变换, 姿态项 */
    Eigen::Matrix3d velanti = Rotation::skewSymmetric(pvacur.vel);
    Eigen::Matrix3d apsi = cnb* velanti;
  
    /* 速度观测值 */
    Eigen::Vector3d tmp9 = cnb * pvacur.vel;
  
    /* 陀螺零偏项 */
    Eigen::Matrix3d ag;
    ag<<0,0,0,0,0,0,0,0,0;
    
      // 构造GNSS位置观测矩阵
      // construct GNSS position measurement matrix
      Eigen::MatrixXd H_gnsspos;
      H_gnsspos.resize(3, Cov_.rows());
      H_gnsspos.setZero();
  
      Eigen::Matrix3d cnb_fill;
      cnb_fill.block(0,0,1,3) = cnb.block(2,0,1,3);
      cnb_fill.block(1,0,1,3) = cnb.block(1,0,1,3);
      cnb_fill.block(2,0,1,3) = cnb.block(0,0,1,3);
      Eigen::Matrix3d apsi_fill;
      apsi_fill.block(0,0,1,3) = -apsi.block(2,0,1,3);
      apsi_fill.block(1,0,1,3) = -apsi.block(1,0,1,3);
      apsi_fill.block(2,0,1,3) = -apsi.block(0,0,1,3);
      Eigen::Matrix3d ag_fill;
      ag_fill.block(0,0,1,3) = -ag.block(2,0,1,3);
      ag_fill.block(1,0,1,3) = -ag.block(1,0,1,3);
      ag_fill.block(2,0,1,3) = -ag.block(0,0,1,3);
     
      H_gnsspos.block(0, V_ID, 3, 3)   = cnb_fill;
      H_gnsspos.block(0, PHI_ID, 3, 3) = apsi_fill;
      H_gnsspos.block(0, BG_ID, 3, 3) = ag;
      // 位置观测噪声阵
      // construct measurement noise matrix
      Eigen::MatrixXd R_gnsspos;
      Eigen::Vector3d odo_std(0.1,0.1,0.1);
      R_gnsspos = odo_std.cwiseProduct(odo_std).asDiagonal();
      if(timestamp_>=641447 && timestamp_<= 641460){
          std::cout<<__FILE__<<","<<__LINE__<<", "<<" H_gnsspos:\n "<<H_gnsspos<<std::endl;
          std::cout<<__FILE__<<","<<__LINE__<<", "<<" R_gnsspos:\n "<<R_gnsspos<<std::endl;
      }
      Eigen::MatrixXd b_pva(3,1);
      b_pva = cnb*pvacur.vel;
      Eigen::MatrixXd dz(3,1);
      dz << b_pva(0)- odo_speed, b_pva(1),b_pva(2);
      // EKF更新协方差和误差状态
      // do EKF update to update covariance and error state
      EKFUpdate(dz, H_gnsspos, R_gnsspos);
      stateFeedback();
    return true;
  }
