# Mahony互补滤波器解析与调优指南

## 目录

1. [算法原理](#算法原理)
2. [代码实现分析](#代码实现分析)
3. [参数调优指南](#参数调优指南)
4. [优化方案](#优化方案)
5. [实际应用经验](#实际应用经验)

---

## 算法原理
参考链接：https://zhuanlan.zhihu.com/p/438724546

### 1. Mahony滤波器基本概念

Mahony互补滤波器是一种用于融合陀螺仪和加速度计数据的姿态估计算法，基于PI控制器设计。

**核心思想：**

- **陀螺仪**：提供高频姿态变化信息，但存在积分漂移
- **加速度计**：提供低频重力方向参考，但有高频噪声
- **融合策略**：高频跟随陀螺仪，低频修正到加速度计

### 2. 数学基础

#### 四元数微分方程

```
dq/dt = 1/2 * q ⊗ ω
```

其中：

- q = [q₀, q₁, q₂, q₃]ᵀ 是四元数
- ω = [0, ωₓ, ωᵧ, ωᵤ]ᵀ 是角速度四元数

#### PI控制器结构

```
ω_corrected = ω_gyro + Kp * error + Ki * ∫error dt
```

---

## 代码实现分析

### 1. 实现步骤详解

#### 步骤1: 参数初始化

```cpp
double sampleFreq = 1.0 / dt;  // 采样频率
// 提取传感器数据
gx = gyroscope(0); gy = gyroscope(1); gz = gyroscope(2);
ax = accelerometer(0); ay = accelerometer(1); az = accelerometer(2);
// 获取当前四元数状态
q0 = navigationSins_->quat(0); q1 = navigationSins_->quat(1);
q2 = navigationSins_->quat(2); q3 = navigationSins_->quat(3);
```

#### 步骤2: 加速度计有效性检查

```cpp
if (!((std::fabs(ax) < 0.000001) && (std::fabs(ay) < 0.000001) && 
      (std::fabs(az) < 0.0000010)) && flag)
```

**作用**: 只有在加速度计数据有效且flag标志位置1时才进行融合

#### 步骤3: 加速度计数据归一化

```cpp
recipNorm = 1.0 / sqrt(ax * ax + ay * ay + az * az);
ax *= recipNorm; ay *= recipNorm; az *= recipNorm;
```

**作用**: 消除重力加速度幅值的影响，只保留方向信息

#### 步骤4: 计算重力向量估计值

```cpp
halfvx = q1 * q3 - q0 * q2;
halfvy = q0 * q1 + q2 * q3;
halfvz = q0 * q0 - 0.5 + q3 * q3;
```

**数学含义**: 四元数旋转矩阵的第三列，表示导航坐标系的z轴（重力方向）在机体坐标系的投影

**坐标系说明**:

`ax, ay, az` 和 `halfvx, halfvy, halfvz` 都在**同一个坐标系——机体坐标系(body frame)**下：

| 变量 | 坐标系 | 物理意义 |
|------|--------|----------|
| ax, ay, az | 机体坐标系 | 加速度计直接测量的比力，归一化后表示**实测重力方向** |
| halfvx, halfvy, halfvz | 机体坐标系 | 根据姿态四元数，将导航系重力向量 [0,0,1]ᵀ 旋转到机体系，表示**估计重力方向** |

两者必须在同一坐标系下，后续的叉积运算才有意义：`error = a_measured × g_estimated`

#### 步骤5: 计算姿态误差

```cpp
halfex = (ay * halfvz - az * halfvy);
halfey = (az * halfvx - ax * halfvz);
halfez = (ax * halfvy - ay * halfvx);
```

**物理意义**: 叉积结果表示需要绕哪个轴旋转多少角度来纠正姿态误差
**数学表达**: error = a_measured × g_estimated

#### 步骤6: 积分反馈（Ki项）

```cpp
if (twoKi > 0.0) {
    integralFbx_ += twoKi * halfex * (1.0 / sampleFreq);
    integralFby_ += twoKi * halfey * (1.0 / sampleFreq);
    integralFbz_ += twoKi * halfez * (1.0 / sampleFreq);
    gx += integralFbx_; gy += integralFby_; gz += integralFbz_;
}
```

**作用**:

- 对姿态误差进行积分，消除陀螺仪的恒定偏差
- 长期补偿陀螺仪零偏，提高姿态估计的长期稳定性

**重要说明**: integralFbx_、integralFby_、integralFbz_ 不是只增大的，因为halfex/y/z可正可负

#### 步骤7: 比例反馈（Kp项）

```cpp
gx += twoKp * halfex;
gy += twoKp * halfey;
gz += twoKp * halfez;
```

**作用**: 对姿态误差进行比例补偿，快速纠正姿态偏差

#### 步骤8: 四元数微分方程积分

```cpp
// 预乘公共因子
gx *= (0.5 * (1.0 / sampleFreq));
gy *= (0.5 * (1.0 / sampleFreq));
gz *= (0.5 * (1.0 / sampleFreq));

// 四元数更新方程（一阶欧拉积分）
qa = q0; qb = q1; qc = q2;  // 保存当前值
q0 += (-qb * gx - qc * gy - q3 * gz);
q1 += (qa * gx + qc * gz - q3 * gy);
q2 += (qa * gy - qb * gz + q3 * gx);
q3 += (qa * gz + qb * gy - qc * gx);
```

**关于0.5系数的说明**:

这里计算的不是纯粹的角度增量，而是**四元数更新所需的半角度增量**：

| 表达式 | 含义 |
|--------|------|
| `1.0 / sampleFreq` | = dt（采样周期） |
| `gx * dt` | 角度增量 Δθₓ（弧度） |
| `gx * 0.5 * dt` | 半角度增量 Δθₓ/2 |

0.5系数来源于四元数微分方程：dq/dt = ½ q ⊗ ω

一阶欧拉积分展开后：q(k+1) = q(k) + (Δt/2) · (q ⊗ ω)

因此预乘0.5是为了匹配四元数更新公式的数学形式。

**数学对应关系**:

| 代码                                     | 数学公式                                    |
| ---------------------------------------- | ------------------------------------------- |
| `q0 += (-qb * gx - qc * gy - q3 * gz)` | q₀ + (Δt/2)(-q₁ωₓ - q₂ωᵧ - q₃ωᵤ) |
| `q1 += (qa * gx + qc * gz - q3 * gy)`  | q₁ + (Δt/2)(q₀ωₓ + q₂ωᵤ - q₃ωᵧ)  |
| `q2 += (qa * gy - qb * gz + q3 * gx)`  | q₂ + (Δt/2)(q₀ωᵧ - q₁ωᵤ + q₃ωₓ)  |
| `q3 += (qa * gz + qb * gy - qc * gx)`  | q₃ + (Δt/2)(q₀ωᵤ + q₁ωᵧ - q₂ωₓ)  |

#### 步骤9: 四元数归一化

```cpp
recipNorm = 1.0 / sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
q0 *= recipNorm; q1 *= recipNorm; q2 *= recipNorm; q3 *= recipNorm;
```

**必要性**: 数值积分会导致四元数长度漂移，必须归一化保持单位长度特性

---

## 参数调优指南

### 1. 当前参数设置

```cpp
twoKp = 1;      // 比例增益
twoKi = 2e-3;   // 积分增益
```

### 2. twoKp (比例增益) 调试经验

| twoKp值 | 特点         | 适用场景               | 缺点     |
| ------- | ------------ | ---------------------- | -------- |
| 0.1-0.5 | 响应慢，稳定 | 高动态运动，有外力干扰 | 收敛慢   |
| 0.5-1.0 | 平衡性能     | 一般车辆导航           | 中等性能 |
| 1.0-2.0 | 快速响应     | 静态或低动态场景       | 易受干扰 |
| >2.0    | 极快响应     | 初始对准               | 很不稳定 |

### 3. twoKi (积分增益) 调试经验

| twoKi值   | 特点     | 适用场景         | 缺点       |
| --------- | -------- | ---------------- | ---------- |
| 1e-4-1e-3 | 慢积分   | 温度变化大的环境 | 零偏估计慢 |
| 1e-3-5e-3 | 标准积分 | 一般应用         | 标准性能   |
| 5e-3-1e-2 | 快积分   | 零偏变化快的场景 | 易过调     |
| >1e-2     | 过快积分 | 不推荐           | 系统不稳定 |

### 4. 实际调试步骤

#### 第一步: 静态测试

```cpp
// 静态时，观察姿态是否稳定收敛
// 目标：5-10秒内收敛到稳定值
twoKp = 1.0;   // 先设中等值
twoKi = 2e-3;  // 先设标准值
```

#### 第二步: 动态测试

- 车辆直线行驶时，观察姿态是否平滑
- 如果振荡：减小twoKp
- 如果漂移：检查twoKi是否合适

#### 第三步: 长期测试

- 观察几小时后是否有长期漂移
- 有漂移：增大twoKi
- 无漂移但响应慢：可以增大twoKp

---

## 优化方案

### 1. 自适应增益调整

#### A. 基于运动状态的自适应

```cpp
if (vehicle_speed < 0.1) {  // 静止时
    twoKp = 2.0;    // 增大比例增益，快速收敛
    twoKi = 5e-3;   // 增大积分增益，估计零偏
} else {  // 运动时
    twoKp = 0.5;    // 减小增益，避免外力干扰
    twoKi = 1e-3;
}
```

#### B. 基于加速度幅值的自适应

```cpp
double acc_magnitude = sqrt(ax*ax + ay*ay + az*az);
double acc_error = fabs(acc_magnitude - 1.0);  // 与重力的偏差

if (acc_error > 0.3) {  // 存在外力加速度
    twoKp *= 0.1;  // 大幅降低增益
    twoKi *= 0.1;
}
```

### 2. 方向性增益

```cpp
// 对不同轴使用不同增益
double twoKpx = twoKp;
double twoKpy = twoKp * 0.8;  // Y轴增益稍小
double twoKpz = twoKp * 1.2;  // Z轴增益稍大

gx += twoKpx * halfex;
gy += twoKpy * halfey; 
gz += twoKpz * halfez;
```

### 3. 积分限幅

```cpp
// 防止积分饱和
const double INTEGRAL_MAX = 0.1;  // 积分上限
integralFbx_ = std::max(-INTEGRAL_MAX, std::min(INTEGRAL_MAX, integralFbx_));
integralFby_ = std::max(-INTEGRAL_MAX, std::min(INTEGRAL_MAX, integralFby_));
integralFbz_ = std::max(-INTEGRAL_MAX, std::min(INTEGRAL_MAX, integralFbz_));
```

### 4. 多模型切换

```cpp
enum MotionState {
    STATIC,      // 静止
    UNIFORM,     // 匀速
    ACCELERATING // 加速
};

MotionState state = detectMotionState();
switch(state) {
    case STATIC:
        twoKp = 2.0; twoKi = 5e-3; break;
    case UNIFORM:
        twoKp = 1.0; twoKi = 2e-3; break;
    case ACCELERATING:
        twoKp = 0.3; twoKi = 1e-3; break;
}
```

---

## 实际应用经验

### 1. 推荐的初始参数

#### 车载应用

```cpp
twoKp = 0.8;    // 适中的响应速度
twoKi = 2e-3;   // 标准的零偏估计速度
```

#### 高精度应用

```cpp
twoKp = 0.5;    // 更稳定
twoKi = 1e-3;   // 更保守的积分
```

### 2. 调试技巧

#### 观察指标

1. **收敛时间**: 静态对准时间（目标5-10秒）
2. **稳态误差**: 长期运行的姿态漂移
3. **动态响应**: 运动过程中的姿态平滑性
4. **抗干扰能力**: 加减速时的姿态稳定性

#### 常见问题及解决

1. **姿态振荡**: twoKp过大，需要减小
2. **收敛慢**: twoKp过小，需要增大
3. **长期漂移**: twoKi不合适，需要调整
4. **加速时发散**: 需要添加外力检测逻辑

### 3. 性能评估方法

#### 静态测试

```cpp
// 记录静态时的姿态标准差
static_roll_std = calculate_std(roll_samples);
static_pitch_std = calculate_std(pitch_samples);
static_yaw_std = calculate_std(yaw_samples);
```

#### 动态测试

```cpp
// 记录运动时的姿态连续性
attitude_jump = max(abs(attitude[i] - attitude[i-1]));
```

---

## 总结

Mahony互补滤波器虽然算法结构相对固定，但有很多优化空间：

1. **参数自适应**: 根据运动状态和环境条件动态调整
2. **多轴独立调整**: 不同轴使用不同参数
3. **积分管理**: 限幅、重置、条件积分
4. **融合其他信息**: 结合轮速、GPS等信息优化参数

关键是要根据具体应用场景和传感器特性进行系统性的测试和调优。
