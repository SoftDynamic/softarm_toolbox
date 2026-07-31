# SoftArm Toolbox

SoftArm Toolbox 是面向连续体机械臂的符号建模、MATLAB 数值代码生成和
Simulink 仿真工具箱。模型由 TOML 文件描述，工具链据此生成 MATLAB 数值函数、
模型清单和 LaTeX 数学手册。

核心动力学方程为

$$
M(q,p)\ddot q+h(q,\dot q,p)=Q(q,\tau_a,w_B,w_e,p),
$$

其中 $q$ 为基座与软臂的联合广义坐标，$p$ 为有序运行时参数，
$\tau_a$ 为软臂广义力，$w_B$ 为机体系基座扳手，$w_e$ 为 NED 世界系末端扳手。

## 1. 工具箱概述

### 1.1 模型族

| 模型 | 单段坐标 | 运动学假设 | 惯性模型 |
| --- | --- | --- | --- |
| PCC | $[b_x,b_y,l]$ | 分段常曲率，可伸缩 | `distributed` 或 `lumped` |
| Euler–Bernoulli Ritz | $[a_x,a_y]$ | 小挠度、小转角、不可伸缩 | 分布惯性 |
| Cosserat PCS | $[\delta\kappa_x,\delta\kappa_y,\delta\kappa_z,\delta\nu_x,\delta\nu_y,\delta\nu_z]$ | 大转角、可剪切、可伸缩、可扭转 | `distributed` 或 `lumped` |

三类模型均支持固定基座和 ZYX 欧拉角浮动基座。浮动基座坐标排列为

```text
[base_x, base_y, base_z, base_roll, base_pitch, base_yaw, arm...]
```

基座刚体、安装变换、软臂和末端负载统一装配到质量矩阵与偏置项中，
从而保留基座—软臂动力学耦合。

### 1.2 功能范围

- PCC 与 Euler–Bernoulli Ritz 多段软臂建模，以及 Cosserat PCS 单段参考模型
- 固定基座与浮动基座联合动力学
- 解析积分与 Gauss–Legendre 积分
- 通用跨段绳索驱动与严格绳长加速度约束
- 平面单点接触、摩擦反力和稳定化
- MATLAB 数值函数、Simulink Model Reference 和 LaTeX 数学手册生成
- 数值线性化与关键节点姿态回放

### 1.3 运行环境

仓库参考模型使用以下版本验证：

```text
Python                3.12.13
MATLAB                R2025a
MATLAB Engine         25.1.2
Wolfram Mathematica   15.0.1
```

Wolfram 后端为可选组件；默认后端为 SymPy。选择 Wolfram 后端时，一个持久
Kernel 会话会覆盖模型推导、解析积分、表达式优化、CSE 和代码生成。模型公式
仍只在 Python/SymPy builder 中维护，所有公共结果仍为 SymPy 表达式。

## 2. 安装与快速开始

### 2.1 安装

```shell
python -m pip install -r requirements.lock
python -m pip install -e .
```

使用 Wolfram 后端时，可将 `.softarm.local.toml.example` 复制为
`.softarm.local.toml`，并填写本机 Wolfram Kernel 路径；也可在构建命令中使用
`--wolfram-kernel` 指定路径。

### 2.2 验证、生成与检查

```shell
softarm validate examples/config/pcc_lumped_n2.toml
softarm build examples/config/pcc_lumped_n2.toml \
  --backend sympy \
  --target matlab \
  --out build/pcc
softarm inspect build/pcc
```

Wolfram 表达式优化后端示例：

```shell
softarm build examples/config/pcc_distributed_n2.toml \
  --backend wolfram \
  --target matlab \
  --out build/pcc-wolfram
```

生成目录包含：

| 文件 | 内容 |
| --- | --- |
| `manifest.json` | 模型族、坐标、参数、执行器和约束元数据 |
| `softarm_*.m` | MATLAB 数值函数 |
| `softarm_model.tex` | 与当前模型对应的 LaTeX 数学手册 |

## 3. 配置参考

模型配置采用 TOML 格式。完整配置位于
[`examples/config`](examples/config/)。

| 配置节 | 用途 |
| --- | --- |
| `[model]` | 模型族、段数和惯性类型 |
| `[base]` | 固定或浮动基座及安装位姿 |
| `[integration]` | 材料坐标积分方法与阶数 |
| `[ritz]` | Euler 模型的归一化 Ritz 多项式 |
| `[parameters]` | 几何、惯性、弹性、阻尼和载荷参数 |
| `[actuation]` | 执行器族、通道及加速度模式 |
| `[constraint]` | 接触或其他加速度级约束 |

### 3.1 模型与积分

PCC 示例：

```toml
[model]
family = "pcc"
segments = 2
inertia = "lumped"

[integration]
method = "analytic"
```

`inertia` 可取 `distributed` 或 `lumped`。数值积分使用
`method = "gauss"` 并设置正整数 `order`。

Cosserat PCS 分布惯性要求显式选择至少二阶 Gauss 积分：

```toml
[model]
family = "cosserat_pcs"
segments = 1
inertia = "distributed"

[integration]
method = "gauss"
order = 2
```

集中惯性使用 `inertia = "lumped"`，并省略 `[integration]`。每段六个坐标是
相对于可配置参考应变的增量，因此零坐标始终表示无应力参考构型。
当前参考包和完整数值验证覆盖单段 Cosserat PCS。配置和 builder 保留通用
`segments` 装配路径，但多段 Cosserat PCS 尚未验证，不作为当前版本的承诺能力。

Euler–Bernoulli Ritz 模型按

$$
\psi(\xi)=\sum_{k=0}^{m}c_k\xi^k
$$

配置两个弯曲平面的模态系数：

```toml
[model]
family = "euler"
segments = 2

[ritz]
x = [0.0, 0.0, 1.5, -0.5]
y = [0.0, 0.0, 1.5, -0.5]
```

系数满足

$$
\psi(0)=\psi'(0)=0,\qquad \psi(1)=1.
$$

### 3.2 基座

固定基座为默认模式。浮动基座配置如下：

```toml
[base]
mode = "floating_rpy"
mount_xyz = [0.0, 0.0, 0.0]
mount_rpy = [0.0, 0.0, 0.0]
```

`mount_xyz` 和 `mount_rpy` 定义基座刚体到软臂根部的固定安装变换。

### 3.3 绳索驱动

绳索通道可跨越任意已声明的软臂段：

```toml
[actuation]
family = "tendon"
acceleration = "strict"

[[actuation.channels]]
name = "t1"
kind = "unilateral"

[[actuation.channels.spans]]
section = 1
radius = 0.0275
angle = 0.0
```

`kind = "unilateral"` 表示非负拉力单绳；
`kind = "signed"` 表示拮抗绳索对的有符号等效通道。

`acceleration = "strict"` 生成绳长加速度约束求解器，适用于通道独立且
通道数不超过软臂坐标数的配置。`acceleration = "none"` 生成力映射，
适用于超驱动或相关通道组。

### 3.4 平面单点接触

```toml
[constraint]
family = "plane_point_contact"
tool_offset = [0.0, 0.0, 0.0]
plane_point = [0.0, 0.0, 0.45]
plane_normal = [0.0, 0.0, -1.0]
friction = 0.3
friction_velocity = 0.01
stabilization_frequency = 20.0
stabilization_ratio = 1.0
```

`plane_normal` 指向自由空间。该约束描述持续接触阶段；法向乘子的可行性输出
可用于上层接触状态管理。

## 4. 生成模型参考

### 4.1 模型清单

`manifest.json` 定义生成模型的公共元数据：

- `model`：模型族、段数和基座模式
- `coordinates.base`、`coordinates.arm`：有序坐标名称
- `parameters`：有序参数名称与默认值
- `actuation`：执行器类型、模式和通道
- `constraint`：约束类型和通道

MATLAB 加载器按清单建立坐标维数、默认参数和可用函数句柄。

### 4.2 MATLAB 数值函数

所有生成包提供以下函数：

```text
softarm_mass(q,p)
softarm_bias(q,dq,p)
softarm_kinematics(q,p)
softarm_end_jacobian(q,p)
softarm_applied_force(q,tau_arm,w_vehicle,w_tip,p)
softarm_forward_dynamics(q,dq,tau_arm,w_vehicle,w_tip,p)
softarm_state_rhs(x,tau_arm,w_vehicle,w_tip,p)
```

配置包含执行器或约束时，生成包还会提供相应的
`softarm_actuator_*` 或 `softarm_constraint_*` 函数。

### 4.3 LaTeX 数学手册

每次构建均生成 `softarm_model.tex`。文档依次给出模型摘要、符号与参数、
运动学、能量与动力学、执行器和约束，并与同一生成包中的数值函数对应。

附加精确符号表达式：

```shell
softarm build examples/config/euler_ritz_n2.toml \
  --target matlab \
  --out build/euler \
  --tex-appendix
```

`--tex-appendix` 将 CAS 优化后的
$V,D,M,h,H_e,J_e,B_v$ 及适用的执行器、约束表达式加入附录。

生成 PDF：

```shell
cd build/euler
pdflatex softarm_model.tex
```

TeX 文档依赖 `article`、`amsmath`、`amssymb`、`geometry` 和 `longtable`。

## 5. MATLAB 接口

```matlab
addpath("matlab")
plant = softarm.loadModel("examples/generated/pcc_lumped_n2");

p = softarm.packParameters(plant.manifest, struct("s1_mass",0.25));
q = [0;0;0.45;0;0;0.50];
dq = zeros(6,1);

ddq = plant.forwardDynamics( ...
    q,dq,zeros(plant.narm,1),zeros(6,1),zeros(6,1),p);

[A,Barm,Bvehicle,Btip] = softarm.linearize( ...
    plant,[q;dq],zeros(plant.narm,1),zeros(6,1),zeros(6,1),p);
```

`softarm.loadModel` 返回的 `plant` 结构包含：

| 字段 | 含义 |
| --- | --- |
| `manifest` | 模型清单 |
| `nbase`、`narm`、`nq` | 基座、软臂和总坐标维数 |
| `parameters` | 默认参数向量 |
| `mass`、`bias`、`kinematics` | 核心模型函数 |
| `endJacobian`、`appliedForce` | 末端 Jacobian 与外力映射 |
| `forwardDynamics`、`stateRhs` | 前向动力学与状态方程 |
| `actuation` | 执行器函数；由配置决定 |
| `constraint` | 约束函数；由配置决定 |

绳索驱动示例：

```matlab
plant = softarm.loadModel( ...
    "examples/generated/pcc_three_tendon_extensible_n2");
T = softarm.packActuatorInputs( ...
    plant.actuation,struct("t1",2,"t2",1.5,"t3",1));
[tau,isFeasible] = plant.actuation.force(q,T,plant.parameters);
ddq = plant.forwardDynamics( ...
    q,dq,tau,zeros(6,1),w,plant.parameters);
```

## 6. Simulink 模型

| 模型 | 用途 |
| --- | --- |
| `softarm_plant.slx` | 无约束 Plant；输入软臂广义力、基座扳手和末端扳手 |
| `softarm_constrained_plant.slx` | 约束 Plant；增加约束加速度并输出反力与可行性 |
| `softarm_actuator_force_block.slx` | 绳索拉力到软臂广义力的 Model Reference |
| `softarm_actuator_acceleration_block.slx` | 严格绳长加速度约束的 Model Reference |
| `softarm_tendon_force_demo.slx` | 绳索拉力驱动示例 |
| `softarm_tendon_acceleration_demo.slx` | 绳长加速度驱动示例 |
| `softarm_flying_contact_demo.slx` | 浮动基座和平面单点接触示例 |

Plant 模型的 `Bundle` System Mask 用于选择生成包，并据此配置基座、软臂、
参数、执行器和约束维数。`backbone_poses` 端口按列堆叠每段末端的世界系
$4\times4$ 齐次变换，输出维数为 $16N\times1$。

在 Demo 的 `Plant` 块中，将 `Bundle directory` 设置为 MATLAB 字符串表达式，
即可切换模型，例如：

```matlab
'examples/generated/euler_two_signed_pairs_n2'
```

### 6.1 姿态回放

三个 Demo 将广义坐标记录到 `softarm_q_log`。仿真结束后运行：

```matlab
softarm_pose_playback
```

也可显式传入日志和生成包：

```matlab
softarm_pose_playback(softarm_q_log)
softarm_pose_playback( ...
    softarm_q_log,Bundle="examples/generated/euler_ritz_n2")
```

播放器显示世界原点、软臂中心线、段末节点和 RGB 姿态轴，并提供时间轴、
播放控制和视角操作。输入可为 timeseries、timetable、structure-with-time、
`Simulink.SimulationOutput` 或首列为时间的数值矩阵。

## 7. 数理基础

### 7.1 符号与坐标约定

第 $i$ 段采用归一化材料坐标 $\xi\in[0,1]$，物理弧长坐标为
$s=L_i\xi$。各段广义坐标按段依次拼接为 $q$，速度为 $\dot q$。齐次变换

$$
{}^0H_i(q,\xi)=
\begin{bmatrix}
{}^0R_i(q,\xi)&{}^0r_i(q,\xi)\\
0&1
\end{bmatrix}
$$

把第 $i$ 段材料截面映射到 NED 基坐标系；多段模型由相邻段变换从基座到
末端依次左乘得到。

材料点的线速度和角速度 Jacobian 定义为

$$
J_{v,i}=\frac{\partial\,{}^0r_i}{\partial q},\qquad
J_{\omega,i}^{(:,j)}=
\operatorname{vex}\!\left(
\operatorname{skew}\!\left(
\frac{\partial\,{}^0R_i}{\partial q_j}\,{}^0R_i^T
\right)\right),
$$

其中 $\operatorname{skew}(A)=(A-A^T)/2$。线性 Euler 模型采用该式的一阶
小转角形式。末端空间 Jacobian 为

$$
J_e(q)=
\begin{bmatrix}
J_{v,e}\\J_{\omega,e}
\end{bmatrix}
\in\mathbb{R}^{6\times n_q}.
$$

### 7.2 PCC 运动学

第 $i$ 段坐标为

$$
q_i=[b_{x,i},\ b_{y,i},\ l_i]^T,\qquad
\theta_i=\sqrt{b_{x,i}^2+b_{y,i}^2}.
$$

$b_x,b_y$ 是弯曲角的笛卡尔分量，直线构型对应
$(b_x,b_y)=(0,0)$。对任意截面 $\xi$，令

$$
x=\xi b_x,\qquad y=\xi b_y,\qquad \rho^2=x^2+y^2,
$$

并定义在原点解析延拓的整函数

$$
S(\rho^2)=\frac{\sin\rho}{\rho},\quad S(0)=1,\qquad
C(\rho^2)=\frac{1-\cos\rho}{\rho^2},\quad C(0)=\frac12.
$$

单段 PCC 变换为

$$
R(\xi)=
\begin{bmatrix}
1-x^2C & -xyC & xS\\
-xyC & 1-y^2C & yS\\
-xS & -yS & 1-\rho^2C
\end{bmatrix},
\qquad
r(\xi)=l\xi
\begin{bmatrix}
xC\\yC\\S
\end{bmatrix}.
$$

该坐标在零曲率处具有唯一表示，并通过 $S$、$C$ 的解析延拓保持运动学及其
导数的连续性。

### 7.3 Euler–Bernoulli Ritz 运动学

每段在两个弯曲平面分别使用一个 Ritz 模态：

$$
w_x(\xi)=a_x\psi_x(\xi),\qquad
w_y(\xi)=a_y\psi_y(\xi),
$$

$$
\psi(\xi)=\sum_{k=0}^{m}c_k\xi^k,\qquad
\psi(0)=\psi'(0)=0,\qquad \psi(1)=1.
$$

局部中心线和截面小转角为

$$
r(\xi)=
\begin{bmatrix}
a_x\psi_x(\xi)\\
a_y\psi_y(\xi)\\
L\xi
\end{bmatrix},
\qquad
\alpha_x=\frac{a_x}{L}\psi_x'(\xi),\quad
\alpha_y=\frac{a_y}{L}\psi_y'(\xi),
$$

$$
R(\xi)\approx
\begin{bmatrix}
1&0&\alpha_x\\
0&1&\alpha_y\\
-\alpha_x&-\alpha_y&1
\end{bmatrix}.
$$

段间变换在拼接后截断到广义坐标的一阶：

$$
H_{\mathrm{lin}}(q)=H(0)+
\sum_j
\left.
\frac{\partial H}{\partial q_j}
\right|_{q=0}q_j.
$$

该模型适用于小挠度、小斜率和小转角工况。Euler 段的弯曲势能为

$$
V_{e,i}=
\frac12\frac{EI_{y,i}a_{x,i}^2}{L_i^3}
\int_0^1\left(\psi_x''\right)^2d\xi
+
\frac12\frac{EI_{x,i}a_{y,i}^2}{L_i^3}
\int_0^1\left(\psi_y''\right)^2d\xi.
$$

### 7.4 Cosserat PCS 运动学

每段使用参考应变 $\kappa_{0,i},\nu_{0,i}$ 和六维应变增量：

$$
\kappa_i=\kappa_{0,i}+\delta\kappa_i,\qquad
\nu_i=\nu_{0,i}+\delta\nu_i.
$$

令 $\Omega_i(\xi)=L_i\xi\widehat{\kappa_i}$、
$z_i(\xi)=(L_i\xi)^2\kappa_i^T\kappa_i$，并定义

$$
\mathcal T(z)=\frac{1-\mathcal S(z)}{z},\qquad \mathcal T(0)=\frac16.
$$

则精确常应变变换为

$$
R_i=I+\mathcal S(z_i)\Omega_i+\mathcal C(z_i)\Omega_i^2,
$$

$$
r_i=\left[I+\mathcal C(z_i)\Omega_i+\mathcal T(z_i)\Omega_i^2\right]L_i\xi\nu_i,
\qquad
H_i(\xi)=\begin{bmatrix}R_i&r_i\\0&1\end{bmatrix}.
$$

参考应变默认为 $\kappa_0=0,\nu_0=[0,0,1]^T$，也可通过运行时参数设置
预弯、预扭或参考剪切。弹性能为

$$
V_{e,i}=\frac{L_i}{2}\delta\xi_i^T
\operatorname{diag}(EI_x,EI_y,GJ,GA_x,GA_y,EA)\delta\xi_i.
$$

### 7.5 能量、质量矩阵与偏置力

对于均匀分布惯性，$m_i$ 和 $I_i$ 表示整段总质量和整段局部惯量。
归一化材料坐标下 $dm=m_i\,d\xi$，每段质量矩阵贡献为

$$
M_i(q)=\int_0^1\left[
m_iJ_{v,i}^TJ_{v,i}
+J_{\omega,i}^T\,{}^0R_iI_i{}^0R_i^T\,J_{\omega,i}
\right]d\xi.
$$

线性 Euler 模型采用参考姿态 ${}^0R_i(0)$ 的惯量旋转，以保持运动学一阶、
能量二阶的一致截断。PCC `lumped` 模式在 $\xi=1/2$ 处计算该段贡献。
末端负载贡献为

$$
M_e=m_eJ_{v,e}^TJ_{v,e}
+J_{\omega,e}^T\,{}^0R_eI_e{}^0R_e^T\,J_{\omega,e}.
$$

总动能为

$$
T(q,\dot q)=\frac12\dot q^TM(q)\dot q.
$$

NED 坐标系中重力沿 $+z$ 方向，重力势能为

$$
V_g=-\sum_i m_i g\int_0^1z_i(q,\xi)d\xi-m_egz_e(q).
$$

PCC `lumped` 模式的对应项为 $-m_i g z_i(q,1/2)$。PCC 的有效弹性势能为

$$
V_{e,i}^{\mathrm{PCC}}=
\frac12k_{bx,i}b_{x,i}^2+
\frac12k_{by,i}b_{y,i}^2+
\frac12k_{l,i}(l_i-L_{0,i})^2.
$$

阻尼采用 Rayleigh 耗散函数

$$
\mathcal{R}(\dot q)=\frac12\dot q^TD\dot q,\qquad
\frac{\partial\mathcal{R}}{\partial\dot q}=D\dot q.
$$

Coriolis 与离心项以能量等价向量表示：

$$
c(q,\dot q)=
\frac{\partial(M\dot q)}{\partial q}\dot q
-\frac12
\left(
\frac{\partial(\dot q^TM\dot q)}{\partial q}
\right)^T.
$$

因此

$$
h(q,\dot q,p)=c(q,\dot q)+\nabla_qV(q,p)+D(p)\dot q,
$$

$$
M(q,p)\ddot q+h(q,\dot q,p)
=Q(q,\tau_a,w_B,w_e,p).
$$

设 $S_a$ 为软臂坐标选择矩阵，$J_B$ 为基座刚体 Jacobian，
$R_{WB}$ 为机体系到 NED 世界系的旋转，则

$$
Q=S_a\tau_a+B_v(q)w_B+J_e^Tw_e,\qquad
B_v=J_B^T\operatorname{diag}(R_{WB},R_{WB}).
$$

外力项由虚功定义：

$$
\delta W=Q^T\delta q.
$$

状态取 $x=[q^T,\dot q^T]^T$，状态方程为

$$
\dot x=
\begin{bmatrix}
\dot q\\
M^{-1}(Q-h)
\end{bmatrix}.
$$

`softarm.linearize` 在指定工作点对该状态方程进行中心差分，返回
$A=\partial\dot x/\partial x$ 以及软臂广义力、基座扳手和末端扳手对应的
输入矩阵。

### 7.6 材料坐标积分

解析积分为默认策略。Gauss–Legendre 积分显式配置为 $n$ 阶时，

$$
\int_0^1f(\xi)d\xi\approx
\frac12\sum_{k=1}^{n}w_k
f\!\left(\frac{\eta_k+1}{2}\right),
$$

其中 $(\eta_k,w_k)$ 为区间 $[-1,1]$ 上的 Gauss–Legendre 节点和权重。
构建过程严格采用配置指定的积分策略。

### 7.7 绳索驱动

第 $a$ 个 PCC 单绳通道的长度坐标为

$$
y_a=\sum_{i\in\mathcal P_a}\left[
l_i-r_{ai}
\left(
b_{x,i}\cos\alpha_{ai}+b_{y,i}\sin\alpha_{ai}
\right)
\right].
$$

`signed` 通道采用相对截面中心对称布置的拮抗绳索半差动长度，
其轴向分量为零。Euler 通道使用 Ritz 端部斜率构造相应弯曲项。

Cosserat PCS 在偏置 $r=[r\cos\theta,r\sin\theta,0]^T$ 处使用精确常应变绳长

$$
\ell_+=L\|\nu+\kappa\times r\|,
$$

单绳通道取 $\ell_+$，signed 对取
$\frac{L}{2}(\|\nu+\kappa\times r\|-\|\nu-\kappa\times r\|)$。

全部通道统一满足

$$
J_a=\frac{\partial y}{\partial q_a},\qquad
\tau_a=-J_a^TT.
$$

严格绳长加速度通过完整系统 Jacobian $J_{a,f}=[0\;J_a]$ 施加：

$$
\begin{bmatrix}
M&J_{a,f}^T\\
J_{a,f}&0
\end{bmatrix}
\begin{bmatrix}
\ddot q\\T
\end{bmatrix}
=
\begin{bmatrix}
Q-h\\
\ddot y_{\mathrm{cmd}}-\dot J_a\dot q_a
\end{bmatrix}.
$$

### 7.8 加速度级约束

约束模型定义约束值、Jacobian、速度偏置、反力映射和稳定化参数：

$$
A(q,p_c)\ddot q+\dot A(q,\dot q,p_c)\dot q=b,
\qquad
Q_c=G(q,\dot q,p_c)\lambda.
$$

理想约束取 $G=A^T$；一般反力模型可定义独立的 $G$。约束动力学采用

$$
\begin{bmatrix}
M&-G\\
A&0
\end{bmatrix}
\begin{bmatrix}
\ddot q\\\lambda
\end{bmatrix}
=
\begin{bmatrix}
Q-h\\b-\dot A\dot q
\end{bmatrix}.
$$

`plane_point_contact` 使用单边法向约束，并以平滑库仑模型将切向摩擦合入
反力映射。求解结果同时提供法向反力可行性和 KKT 条件诊断。

## 8. 示例索引

| 示例 | 配置与生成包 |
| --- | --- |
| 两段 PCC，中点集中惯性 | `pcc_lumped_n2` |
| 两段 PCC，分布惯性 | `pcc_distributed_n2` |
| 两段 Euler–Bernoulli Ritz | `euler_ritz_n2` |
| 两段可伸缩 PCC，三绳驱动 | `pcc_three_tendon_extensible_n2` |
| 两段 PCC，signed 绳索对 | `pcc_signed_pair_n2` |
| 两段 Euler，两个 signed 绳索对 | `euler_two_signed_pairs_n2` |
| 浮动基座 PCC，平面单点接触 | `pcc_flying_plane_contact_n1` |
| 单段 Cosserat PCS，分布惯性与三绳驱动 | `cosserat_pcs_distributed_tendon_n1` |
| 单段 Cosserat PCS，中点集中惯性 | `cosserat_pcs_lumped_n1` |

MATLAB 绳索驱动脚本见
[`examples/actuation`](examples/actuation/)。
