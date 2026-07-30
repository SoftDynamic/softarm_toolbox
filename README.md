# SoftArm Toolbox

这是一个以 SymPy 为唯一数学表达式来源的连续体机械臂建模、MATLAB 数值仿真和控制设计基础工具箱。Wolfram 15 是可选重型 CAS；MATLAB/Simulink 只消费生成的数值函数，不参与符号建模。

## 模型

- **PCC**：每段坐标为 `[bx, by, l]`。`bx/by` 是正交弯曲角分量，避免 `[theta,phi]` 在直线构型下的方位角奇异。支持沿材料坐标积分的 `distributed` 惯性和中点质量 `lumped` 惯性。
- **Euler–Bernoulli Ritz**：每段坐标为 `[ax, ay]`，采用小挠度线性假设模态。配置必须显式给出两个归一化 Ritz 多项式；该模型不适合段内大斜率。

两种软臂模型都可使用固定基或 ZYX 欧拉角浮动基。浮动基坐标
`[base_x,base_y,base_z,base_roll,base_pitch,base_yaw]` 位于软臂坐标之前，无人机刚体、安装变换、软臂和末端负载由同一能量表达式装配，因此质量矩阵包含完整的基座—软臂耦合项。联合方程为

$$
M(q,p) \ddot{q} + h(q,\dot{q},p) = Q(q,\tau_a,w_B,w_e,p)
$$

$\tau_a$ 是仅对应软臂自由度的抽象广义力；$w_B$ 是作用在无人机质心、在机体系表达的扳手；$w_e$ 是在 NED 世界系表达的末端扳手。`softarm_applied_force` 将三者统一映射到完整广义力。实际绳索或气动执行器通过独立的 `ActuatorMap` 映射到 $\tau_a$。

## 数理基础

### 符号与坐标约定

第 $i$ 段采用归一化材料坐标 $\xi\in[0,1]$，物理弧长坐标为
$s=L_i\xi$。各段广义坐标按段依次拼接为 $q$，速度为
$\dot q$。齐次变换

$$
{}^0H_i(q,\xi)=
\begin{bmatrix}
{}^0R_i(q,\xi)&{}^0r_i(q,\xi)\\
0&1
\end{bmatrix}
$$

把第 $i$ 段材料截面映射到 NED 基坐标系；多段模型由相邻段变换从基座到末端依次左乘得到。

材料点的线速度和角速度 Jacobian 定义为

$$
J_{v,i}=\frac{\partial\,{}^0r_i}{\partial q},\qquad
J_{\omega,i}^{(:,j)}=
\operatorname{vex}\!\left(
\operatorname{skew}\!\left(
\frac{\partial\,{}^0R_i}{\partial q_j}\,{}^0R_i^T
\right)\right),
$$

其中 $\operatorname{skew}(A)=(A-A^T)/2$。线性 Euler 模型使用该式的一阶小转角形式。末端空间 Jacobian 为

$$
J_e(q)=\begin{bmatrix}J_{v,e}\\J_{\omega,e}\end{bmatrix}
\in\mathbb{R}^{6\times n_q}.
$$

### PCC 无奇异弯曲坐标与运动学

第 $i$ 段坐标为

$$
q_i=[b_{x,i},\ b_{y,i},\ l_i]^T,\qquad
\theta_i=\sqrt{b_{x,i}^2+b_{y,i}^2}.
$$

$b_x,b_y$ 是弯曲角的笛卡尔分量，因此直线构型由
$(b_x,b_y)=(0,0)$ 唯一表示，不需要不可观的弯曲方位角。对任意截面 $\xi$，令

$$
x=\xi b_x,\qquad y=\xi b_y,\qquad \rho^2=x^2+y^2,
$$

并定义在原点解析延拓的整函数

$$
S(\rho^2)=\frac{\sin\rho}{\rho},\quad S(0)=1,\qquad
C(\rho^2)=\frac{1-\cos\rho}{\rho^2},\quad C(0)=\frac12.
$$

单段 PCC 变换使用

$$
R(\xi)=
\begin{bmatrix}
1-x^2C & -xyC & xS\\
-xyC & 1-y^2C & yS\\
-xS & -yS & 1-\rho^2C
\end{bmatrix},
\qquad
r(\xi)=l\xi
\begin{bmatrix}xC\\yC\\S\end{bmatrix}.
$$

`SincSqrt` 和 `CoscSqrt` 在 SymPy 中直接表示上述两个整函数，并提供解析的一、二阶导数。生成的 MATLAB 函数只在
$|\rho^2|<10^{-8}$ 时切换到 Taylor 多项式；这不会修改 $q$，也不会在质量矩阵中加入人为正则项。

### Euler–Bernoulli Ritz 模型

Euler 模型每段每个弯曲平面使用一个显式配置的 Ritz 模态：

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

段间变换在拼接后严格截断到广义坐标的一阶：

$$
H_{\mathrm{lin}}(q)=H(0)+
\sum_j\left.\frac{\partial H}{\partial q_j}\right|_{q=0}q_j.
$$

因此该模型是线性小挠度梁模型，不应被用于段内大斜率或大转角。若需要这类行为，应增加几何精确梁或
Cosserat 模型，而不是在当前 Ritz 模型中混入部分高阶项。

Euler 段的弯曲势能为

$$
V_{e,i}=\frac12\frac{EI_{y,i}a_{x,i}^2}{L_i^3}
\int_0^1\left(\psi_x''\right)^2d\xi
+\frac12\frac{EI_{x,i}a_{y,i}^2}{L_i^3}
\int_0^1\left(\psi_y''\right)^2d\xi.
$$

### 能量、质量矩阵与偏置力

对于均匀分布惯性，配置中的 $m_i$ 和 $I_i$ 表示整段总质量和整段局部惯量；由于使用归一化材料坐标，
$dm=m_i\,d\xi$。每段质量矩阵贡献为

$$
M_i(q)=\int_0^1\left[
m_iJ_{v,i}^TJ_{v,i}
+J_{\omega,i}^T\,{}^0R_iI_i{}^0R_i^T\,J_{\omega,i}
\right]d\xi.
$$

对线性 Euler 模型，为保持运动学一阶、能量二阶的一致截断，上式的惯量旋转取参考姿态
${}^0R_i(0)$。PCC `lumped` 模式用 $\xi=1/2$ 处的同一表达式替代积分。末端负载贡献为

$$
M_e=m_eJ_{v,e}^TJ_{v,e}
+J_{\omega,e}^T\,{}^0R_eI_e{}^0R_e^T\,J_{\omega,e}.
$$

总动能为

$$
T(q,\dot q)=\frac12\dot q^TM(q)\dot q.
$$

NED 坐标系中重力沿 $+z$ 方向，故重力势能写为

$$
V_g=-\sum_i m_i g\int_0^1 z_i(q,\xi)d\xi-m_egz_e(q).
$$

PCC `lumped` 模式相应地以 $-m_i g z_i(q,1/2)$ 代替该段的重力积分。

PCC 的有效弹性势能为

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

代码不显式构造非唯一的 Coriolis 矩阵 $C$，而是直接生成能量等价向量

$$
c(q,\dot q)=
\frac{\partial(M\dot q)}{\partial q}\dot q
-\frac12
\left(\frac{\partial(\dot q^TM\dot q)}{\partial q}\right)^T.
$$

于是

$$
h(q,\dot q,p)=c(q,\dot q)+\nabla_qV(q,p)+D(p)\dot q,
$$

$$
M(q,p)\ddot q+h(q,\dot q,p)=Q(q,\tau_a,w_B,w_e,p).
$$

外力项来自虚功

$$
\delta W=Q^T\delta q.
$$

状态取 $x=[q^T,\dot q^T]^T$，`softarm_state_rhs` 输出

$$
\dot x=
\begin{bmatrix}
\dot q\\
M^{-1}(Q-h)
\end{bmatrix}.
$$

MATLAB 的 `softarm.linearize` 在给定工作点对该状态方程做中心差分，返回
$A=\partial\dot x/\partial x$ 以及软臂、无人机扳手和末端扳手三组输入矩阵。

### 解析积分与显式数值积分

默认由所选 CAS 求解材料坐标上的解析积分。若配置显式选择 $n$ 阶 Gauss–Legendre 积分，则

$$
\int_0^1f(\xi)d\xi\approx
\frac12\sum_{k=1}^{n}w_k
f\!\left(\frac{\eta_k+1}{2}\right),
$$

其中 $(\eta_k,w_k)$ 是区间 $[-1,1]$ 上的 Gauss–Legendre 节点和权重。解析失败不会自动降级为数值积分。

## 环境与安装

本仓库已在以下环境验证：

```text
Python  3.12.13
MATLAB  R2025a / matlabengine 25.1.2
Wolfram Mathematica 15.0.1
```

安装已锁定依赖：

```shell
python -m pip install -r requirements.lock
python -m pip install -e .
```

复制 `.softarm.local.toml.example` 为 `.softarm.local.toml` 可覆盖本机工具路径；该文件不会提交。

## 配置与生成

参考配置位于 `examples/config/`。Euler Ritz 系数按
`psi(xi)=sum(c[k]*xi^k)` 给出，必须满足 `c0=c1=0` 和 `sum(c)=1`。

```shell
softarm validate examples/config/pcc_lumped_n2.toml
softarm build examples/config/pcc_lumped_n2.toml --backend sympy --target matlab --out build/pcc
softarm build examples/config/pcc_distributed_n2.toml --backend wolfram --target matlab --out build/pcc-wolfram
softarm inspect build/pcc
```

### LaTeX 数学文档

每次 `softarm build` 都会在输出 bundle 中生成固定文件名
`softarm_model.tex`。它是独立的英文 `article` 文档，与 MATLAB 代码共用同一个
`SymbolicPlant`、`ActuationModel` 和 `ConstraintModel`，因此不是另一套手工维护的公式。
正文按符号与参数、运动学、能量与动力学、执行器和约束分层组织；大型质量矩阵和偏置向量
默认只以装配公式表示，适合人工核对并作为论文推导的起点。

需要查看 CAS 生成的精确表达式时，使用：

```shell
softarm build examples/config/euler_ritz_n2.toml --target matlab --out build/euler --tex-appendix
```

`--tex-appendix` 会在同一文件末尾按输出块加入优化和 CSE 后的
$V,D,M,h,H_e,J_e,B_v$，以及适用的执行器和约束表达式。大型浮动基模型的附录可能很长，
普通构建和仓库内参考 bundle 因而不启用此选项。论文中建议引用正文的命名装配公式，并仅把
精确附录作为模型复现或人工核验材料。

本地生成 PDF：

```shell
cd build/euler
pdflatex softarm_model.tex
```

该 TeX 仅依赖标准 `article`、`amsmath`、`amssymb`、`geometry` 和 `longtable`，
不被 MATLAB loader 或 Simulink 消费，且不会改变 manifest。

浮动基和平面接触可直接在同一 TOML 中声明；完整示例见
`examples/config/pcc_flying_plane_contact_n1.toml`：

```toml
[base]
mode = "floating_rpy"
mount_xyz = [0.0, 0.0, 0.0]
mount_rpy = [0.0, 0.0, 0.0]

[parameters]
vehicle_mass = 1.5
vehicle_Ixx = 0.03
vehicle_Iyy = 0.03
vehicle_Izz = 0.05

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

平面法向指向自由空间。该接触实现假定约束已经激活，适用于持续贴合阶段；接近、冲击和脱离应由上层模式管理器根据间隙和反力可行性切换。

解析积分是默认策略。若 CAS 无法完成，构建会失败；只有配置明确指定
`method="gauss"` 和积分阶数时才会生成 Gauss–Legendre 求和表达式。

生成包的 `manifest.json` 是唯一元数据源，包含 `model`、分组的 `coordinates.base/arm`、有序参数、执行器和约束描述。它不包含版本字段，也不兼容旧 manifest。固定 MATLAB 接口为：

```text
softarm_mass(q,p)
softarm_bias(q,dq,p)
softarm_kinematics(q,p)
softarm_end_jacobian(q,p)
softarm_applied_force(q,tau_arm,w_vehicle,w_tip,p)
softarm_forward_dynamics(q,dq,tau_arm,w_vehicle,w_tip,p)
softarm_state_rhs(x,tau_arm,w_vehicle,w_tip,p)
```

## 更新与增加数理模型

### 修改现有模型

1. 在 `src/softarm/geometry.py` 修改单段运动学，所有公式必须由 SymPy 对象构造；不要在 MATLAB 或 Wolfram 文件中复制一套公式。
2. 在 `src/softarm/derive.py` 修改参数、能量或惯性装配。保持 `SymbolicPlant` 的公共结果为 `mass`、`potential`、`damping`、`kinematics` 和 `end_jacobian`；`bias` 会由统一能量公式生成。
3. 在 `src/softarm/config.py` 增加必要的配置校验，并同步更新 `examples/config/`。结构参数进入 TOML；仿真中需要调节的物理参数应进入有序运行时参数向量 `p`。
4. 若修改生成函数的固定签名或输出维度，同步更新 `src/softarm/codegen/matlab.py`、MATLAB 加载层和 `matlab/build_softarm_plant.m`；否则不需要手工编辑 SLX。
5. 重新生成受影响的 `examples/generated/` 参考包，并运行 Python、MATLAB 和 Simulink 测试。

### 增加新的几何/梁模型

1. 编写一个接收 `ModelConfig`、返回 `SymbolicPlant` 的 builder。仓库内置模型加入 `derive.py` 的 `_MODEL_BUILDERS`；外部 Python 代码在调用 `derive` 前执行 `softarm.register_model(name, builder)`。当前 CLI 不自动发现外部插件，如需直接从 CLI 使用，还必须在 CLI 启动阶段显式导入注册模块。
2. 优先复用 `derive.py` 中的通用能量装配思路；若新模型需要新的结构配置字段，应同时扩展 `ModelConfig` 和 TOML 解析，不能把结构常量隐藏在生成器中。
3. 对模型的适用范围作明确选择，例如小挠度/大转角、可伸长/不可伸长、有剪切/无剪切；不要把不同阶次或互相矛盾的假设混入同一模型。
4. 至少增加以下验证：零构型或参考构型、运动学有限差分、质量矩阵对称与名义正定、独立能量/积分计算、Coriolis 能量恒等式、外力虚功，以及 N=1/N=2 多段组合。

### 增加特殊函数或 CAS 节点

若新模型出现可去奇点或 SymPy/Wolfram/MATLAB 不能共同表示的函数：

1. 在 `src/softarm/special.py` 定义 SymPy 函数、所需阶数的解析导数和独立数值实现。
2. 在 `src/softarm/ast.py`、`src/softarm/backends/wolfram_bridge.wls` 和 MATLAB printer 中增加同一节点的往返与输出规则。
3. 为奇点值、邻域连续性、导数和 SymPy/Wolfram/MATLAB 数值一致性增加测试；不得用修改状态变量的 epsilon 技巧代替解析延拓。

### 增加执行器或约束

基础 Plant 分别接收软臂广义力、无人机机体系扳手和末端世界系扳手。绳索、气腔或电机模型应实现 `ActuatorMap`，在 Plant 外部计算

$$
\tau_a=B(q_a,p_u)u
$$

加速度级约束通过 `ConstraintModel` 描述约束值、Jacobian、速度偏置、反力映射和稳定化参数：

$$
A(q,p_c)\ddot q+\dot A(q,\dot q,p_c)\dot q=b,\qquad Q_c=G(q,\dot q,p_c)\lambda.
$$

内置求解器使用无正则化 KKT 方程；理想约束取 $G=A^T$，非理想反力可提供独立的 $G$。`plane_point_contact` 示例生成单边法向约束，并用平滑库仑模型把滑动摩擦合入反力映射。求解器只处理已激活约束；负法向乘子通过可行性输出报告，不自动执行碰撞、接触或脱离切换。

## 通用绳索驱动

可在模型 TOML 中声明任意数量、任意跨段的绳索通道。`unilateral` 表示拉力必须非负的真实单绳；`signed` 表示一对拮抗绳索的理想差动通道，其等效拉力可正可负且不产生轴向作用。

```toml
[actuation]
family = "tendon"
acceleration = "strict" # 或 "none"

[[actuation.channels]]
name = "t1"
kind = "unilateral"
[[actuation.channels.spans]]
section = 1
radius = 0.0275
angle = 0.0
```

第 $a$ 个 PCC 单绳通道的长度坐标为

$$
y_a=\sum_{i\in\mathcal P_a}\left[
l_i-r_{ai}(b_{x,i}\cos\alpha_{ai}+b_{y,i}\sin\alpha_{ai})
\right].
$$

`signed` 通道去掉 $l_i$ 项，等价于相对截面中心对称布置的一对绳索的半差动长度。Euler 使用同一走线配置，但将弯曲项替换为 Ritz 端部斜率；Euler 本身仍不可轴向伸缩。

所有通道统一使用

$$
J_a=\frac{\partial y}{\partial q},\qquad
\tau_a=-J_a^\top T.
$$

`acceleration="strict"` 仅允许独立且不多于软臂坐标数的通道，并生成无正则化 KKT 求解器；`acceleration="none"` 只生成力映射，可用于超驱动或相关绳索组。半径作为有序运行时参数写入现有 `p`，通道名称与类型直接来自 manifest。

```matlab
plant = softarm.loadModel("examples/generated/pcc_three_tendon_extensible_n2");
T = softarm.packActuatorInputs(plant.actuation,struct("t1",2,"t2",1.5,"t3",1));
[tau,isFeasible] = plant.actuation.force(q,T,plant.parameters);
ddq = plant.forwardDynamics(q,dq,tau,zeros(6,1),w,plant.parameters);
```

三绳可伸缩、严格绳长加速度和 signed 等效配对示例见 `examples/actuation/README.md`。特殊滑轮或非标准长度公式通过 `softarm.register_actuator(name,builder)` 注册；builder 必须返回完全由 SymPy 表达式构成的 `ActuationModel`。

## MATLAB 与 Simulink

```matlab
addpath("matlab")
plant = softarm.loadModel("examples/generated/pcc_lumped_n2");
p = softarm.packParameters(plant.manifest, struct("s1_mass",0.25));
q = [0;0;0.45;0;0;0.50];
dq = zeros(6,1);
ddq = plant.forwardDynamics(q,dq,zeros(plant.narm,1),zeros(6,1),zeros(6,1),p);
[A,Barm,Bvehicle,Btip] = softarm.linearize(plant,[q;dq],zeros(plant.narm,1),zeros(6,1),zeros(6,1),p);
```

`softarm_plant.slx` 是无约束底层 Plant，输入为 `tau_arm`、`vehicle_wrench` 和 `tip_wrench`。`softarm_constrained_plant.slx` 增加约束加速度输入，并输出反力、可行性和 KKT 条件诊断。两个模型都带 `Bundle` System Mask，切换 bundle 会自动更新基座、软臂、执行器、约束和参数维数。

模型可由 `matlab/build_softarm_plant.m` 完整重建，不需要手工编辑 SLX。

### 图形化绳索驱动

以下模型可直接从文件浏览器打开，无需预先运行初始化脚本；回调会自动加载 `matlab/` 和默认的三绳可伸缩 PCC 包。

- `softarm_actuator_force_block.slx`：可复用拉力到广义力 Model Reference，输入 `q/tension`，输出 `tau/is_feasible`。
- `softarm_actuator_acceleration_block.slx`：可复用严格绳长加速度 Model Reference，输出总 `tau`、约束拉力、可行性、约束条件数和 `ddq`。
- `softarm_tendon_force_demo.slx`：已经把拉力驱动器与 Plant 串联；可设置拉力、无人机扳手和末端扳手。
- `softarm_tendon_acceleration_demo.slx`：已经把严格加速度驱动器与 Plant 串联；可设置绳长加速度、附加软臂力和两类扳手。
- `softarm_flying_contact_demo.slx`：浮动基 PCC 与平面单点摩擦约束示例；从已贴合状态开始。

两个 Demo 的 `State scope` 显示 `q/dq`。反馈上的 Memory 块用于跨 Model Reference 打断 Simulink 的保守代数环判断。构建脚本为 `matlab/build_softarm_actuator_models.m`。

切换模型时，双击 Demo 中的 `Plant` 块，在 Mask 的 `Bundle directory` 中填写生成包路径并点击 Apply，例如：

```matlab
'examples/generated/euler_two_signed_pairs_n2'
```

这会在不更改连线的情况下，把三绳可伸缩 PCC 切换为 N=2 Euler + 两个 signed tendon pair；`q/dq` 从 6 维变为 4 维，tendon 输入从 3 维变为 2 维。该字段是 Simulink 模型参数表达式，因此字符串路径需要单引号。自建顶层模型时，引用 `softarm_plant` 的 Model block 会继承相同的 Mask。

## 测试

```shell
python -m pytest
matlab -batch "addpath('matlab'); r=runtests('matlab/tests'); assertSuccess(r)"
```

测试覆盖配置约束、PCC 零曲率的一二阶导数、浮动基质量耦合、Ritz 边界条件、通用绳索虚功/轴向模态、平面约束与摩擦耗散、新 manifest、MATLAB 参考包、数值线性化及 Simulink 短时运行。
