# SoftArm Toolbox

SoftArm Toolbox 是面向连续体机械臂的符号建模、MATLAB/CasADi 数值代码生成和 Simulink 仿真工具箱。模型由 TOML 文件描述，工具链据此生成数值函数、模型清单 和 LaTeX 数学手册。

核心动力学方程为

$$
M(q,p)\ddot q+h(q,\dot q,p)=Q(q,\tau_a,w_B,w_e,p),
$$

其中 $q$ 为基座与软臂的联合广义坐标，$p$ 为有序运行时参数， $\tau_a$ 为软臂广义力，$w_B$ 为机体系基座扳手，$w_e$ 为 NED 世界系末端扳手。

## 1. 工具箱概述

### 1.1 梁模型与空间参数化

| 梁模型 | Ritz | PCS | [PAC](https://arxiv.org/abs/2211.10188) | 保留的应变分量 |
| --- | --- | --- | --- | --- |
| Euler–Bernoulli | ✅ | ✅ | ✅ | $\boldsymbol e=[\kappa_x,\kappa_y,0,0,0,0]^T$ |
| Extensible Euler–Bernoulli | ✅ | ✅ | ✅ | $\boldsymbol e=[\kappa_x,\kappa_y,0,0,0,\varepsilon_z]^T$ |
| Cosserat | ❌ | ✅ | ❌ | $\boldsymbol e=[\kappa_x,\kappa_y,\kappa_z,\gamma_x,\gamma_y,\varepsilon_z]^T$ |

表中

$$
\boldsymbol e=
[\kappa_x,\kappa_y,\kappa_z,\gamma_x,\gamma_y,\varepsilon_z]^T
$$

表示相对于参考构型的杆应变：$\kappa$ 为角应变，$\gamma_x,\gamma_y$ 为横向剪切应变， $\varepsilon_z$ 为轴向应变。该列表示各梁理论保留的应变分量，不表示所有空间参数化都精确 实现非线性几何约束。PCS 和 PAC 精确实现 Euler–Bernoulli 的无剪切约束；不可伸长模型 固定轴向应变，可伸长模型则显式表示轴向伸长。Ritz 使用一阶小挠度 运动学。PAC 的段级截面转角 $\phi_i$ 不是连续扭转应变 $\kappa_z$。

三类模型均支持固定基座和 ZYX 欧拉角浮动基座。浮动基座坐标排列为

```text
[base_x, base_y, base_z, base_roll, base_pitch, base_yaw, arm...]
```

基座刚体、安装变换、软臂和末端负载统一装配到质量矩阵与偏置项中， 从而保留基座—软臂动力学耦合。

上述基座能力同时适用于 `symbolic_lagrange` 和 `recursive`。PCS、PAC 与 Cosserat PCS 使用几何精确的逐段递推；一阶 displacement-Ritz 使用专用 affine jet，在浮动基座下保持基座姿态非线性和基座—软臂耦合，同时严格丢弃软臂坐标间 的二阶项。工具箱尚未实现 strain-Ritz。

### 1.2 功能范围

- Euler–Bernoulli、Extensible Euler–Bernoulli 与 Cosserat rod 建模
- Ritz、分段恒应变（PCS）与三维分段仿射曲率（PAC）空间参数化
- 固定基座与浮动基座联合动力学
- 解析积分与 Gauss–Legendre 积分
- 通用跨段绳索驱动与严格绳长加速度约束
- `symbolic_lagrange` 与 `recursive` 下的平面单点接触、摩擦反力和稳定化
- MATLAB 数值函数、CasADi bundle、Simulink Model Reference 和 LaTeX 数学手册生成
- 数值线性化与关键节点姿态回放

### 1.3 运行环境

仓库参考模型使用以下版本验证：

```text
Python                3.12.13
MATLAB                R2025a
CasADi                3.7.2
Wolfram Mathematica   15.0.1
```

Wolfram 后端为可选组件；默认后端为 SymPy。模型公式、几何装配、积分、执行器 和约束始终由 Python/SymPy 构造。选择 Wolfram 后端时，一个持久 Kernel 负责 bias 的批量求导；`FactorTerms` 和实验性 Wolfram CSE 必须分别通过命令行开关 显式启用。所有公共结果仍为 SymPy 表达式，任何已选择的 Wolfram 操作失败都会 终止构建，不会静默切换算法。

## 2. 安装与快速开始

### 2.1 安装

```shell
uv sync --locked
```

uv 根据 `.python-version` 使用最新可用的 Python 3.12 补丁版本；开发依赖组默认同步。 可将 `.softarm.local.toml.example` 复制为 `.softarm.local.toml`，填写本机 MATLAB、Wolfram Kernel 和 `pdflatex` 路径。Wolfram Kernel 也可在构建命令中使用 `--wolfram-kernel` 指定；仓库维护命令和 TeX 编译测试优先使用本机配置，避免误用 系统中同名但版本不同的工具。

### 2.2 验证、生成与检查

```shell
uv run --locked softarm validate \
  examples/config/extensible_euler_bernoulli_pcs_lumped_n2.toml
uv run --locked softarm build \
  examples/config/extensible_euler_bernoulli_pcs_lumped_n2.toml \
  --backend sympy \
  --target matlab \
  --out build/extensible-euler-bernoulli-pcs
uv run --locked softarm inspect build/extensible-euler-bernoulli-pcs
```

使用 Wolfram 加速 bias 批量求导：

```shell
uv run --locked --extra wolfram softarm build \
  examples/config/extensible_euler_bernoulli_pcs_distributed_n2.toml \
  --backend wolfram \
  --target matlab \
  --out build/extensible-euler-bernoulli-pcs-wolfram
```

成功的 `softarm build` 默认不输出内容。增加 `--verbose` 可输出适合直接复制到 问题报告中的完整构建诊断：

```shell
uv run --locked --extra wolfram softarm build \
  examples/config/extensible_euler_bernoulli_pcs_distributed_n2.toml \
  --backend wolfram \
  --target matlab \
  --out build/extensible-euler-bernoulli-pcs-wolfram \
  --verbose
```

```text
SoftArm Build Diagnostics
=========================

Configuration
-------------
...

Symbolic Strategy
-----------------
...

Build Stages
------------
...

Resolved Model
--------------
...

Exported Symbolic Functions
---------------------------
...

Generated Artifacts
-------------------
...

Build Result
------------
SUCCESS
```

报告使用解析后的绝对配置/输出路径和运行时实际可用的 SymPy、Wolfram 信息； 无法可靠取得的版本或路径不会显示。配置按模型、基座、积分、Ritz、执行器和 约束分组，并列出推导后包含默认值的完整运行时参数。

核心计时从 `build_bundle` 开始，不包括配置文件读取。`MATLAB generation/CSE` 汇总所选 `FactorTerms`、CSE、MATLAB/TeX/manifest 生成。随后执行的诊断分析不计入 `Total build`，并单独报告耗时。每个符号导出函数会显示输出形状、CSE 前后 `sympy.count_ops`、临时量数量和 MATLAB 文件大小；这些运算数是符号复杂度估计， 不是硬件 FLOP。若构建失败，已完成阶段保留，失败阶段标记 `(failed)`，错误原因 仍写到标准错误流。

### 2.3 Wolfram 优化开关

默认情况下，代码生成使用稳定的 `sympy.cse`，且不运行 `FactorTerms`。以下 开关会改变符号执行路径，因此只能与 `--backend wolfram` 一起显式使用：

| 开关 | 作用 | 建议使用场景 |
| --- | --- | --- |
| `--wolfram-factor-terms` | 每个 MATLAB 函数在 CSE 前执行 Wolfram `FactorTerms` | 生成表达式明显过度展开时单独试用；可能显著变慢，也不保证减少临时量 |
| `--wolfram-cse` | 用实验性 Experimental OptimizeExpression 代替 `sympy.cse` | SymPy CSE 已确认是主要瓶颈，或需要比较生成代码时单独试用 |
| `--wolfram-timeout SECONDS` | 设置每项 Wolfram 操作的超时，默认 600 秒 | 已确认操作合理但默认时间不足时调整 |

首次使用建议只选择 `--backend wolfram`。若仍需优化生成代码，分别试用两个 开关并比较构建时间和 MATLAB 文件大小，不建议一开始同时开启。开关启用后， 超时、不支持的实验性返回结构、非法临时量依赖或等价性校验失败都会直接报错。 错误信息会给出关闭相应开关或调整超时的建议，但工具不会自动重试。

例如只试用 `FactorTerms`：

```shell
uv run --locked --extra wolfram softarm build \
  examples/config/extensible_euler_bernoulli_pcs_lumped_n2.toml \
  --backend wolfram \
  --wolfram-factor-terms \
  --out build/extensible-euler-bernoulli-pcs-factor-terms
```

例如只试用实验性 Wolfram CSE：

```shell
uv run --locked --extra wolfram softarm build \
  examples/config/euler_bernoulli_ritz_n2.toml \
  --backend wolfram \
  --wolfram-cse \
  --out build/euler-wolfram-cse
```

参考性能如下：

| 配置 | SymPy 模型推导 | Kernel 启动 | Wolfram bias | MATLAB 生成/CSE | 完整构建 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Euler N2 | 0.347 s | 7.302 s | 0.125 s | 0.067 s | 7.848 s |
| Extensible Euler–Bernoulli PCS lumped N2 | 0.869 s | 7.138 s | 21.977 s | 2.810 s | 32.970 s |
| Extensible Euler–Bernoulli PCS distributed N2 | 1.556 s | 7.407 s | 58.440 s | 5.804 s | 73.256 s |

不同机器上的 Kernel 启动、许可证检查和 CSE 时间会变化，应在自己的构建环境 分别测量。尤其不要因为某个配置较快，就同时打开两个实验开关。

生成目录包含：

| 文件 | 内容 |
| --- | --- |
| `manifest.json` | Rod、参数化、坐标、参数、执行器和约束元数据 |
| `softarm_*.m` | MATLAB 数值函数 |
| `softarm_model.tex` | 与当前模型对应的 LaTeX 数学手册 |

## 3. 配置参考

模型配置采用 TOML 格式。完整配置位于 [`examples/config`](examples/config/)。

| 配置节 | 用途 |
| --- | --- |
| `[model]` | Rod 理论、空间参数化、段数和惯性类型 |
| `[dynamics]` | 必填的动力学装配方法：`symbolic_lagrange` 或 `recursive` |
| `[base]` | 固定或浮动基座及安装位姿 |
| `[integration]` | 材料坐标积分方法与阶数 |
| `[ritz]` | Ritz 参数化的归一化多项式 |
| `[parameters]` | 几何、惯性、弹性、阻尼和载荷参数 |
| `[actuation]` | 执行器族、通道及加速度模式 |
| `[constraint]` | 约束（实验性功能） |

每个配置都必须显式选择动力学装配方法，不提供默认值。`symbolic_lagrange` 保留全局符号欧拉–拉格朗日推导；`recursive` 生成逐段局部核，并在 MATLAB 运行时执行前向运动传播和反向力旋量回拉。递推模式不接受 Wolfram 后端选项或 符号 TeX appendix。内置平面点接触由递推工具点运动学 jet 构造，因此 MATLAB 和 CasADi 均提供与符号模型一致的约束函数和 active-contact KKT。现有 displacement-Ritz 由保持全局一阶小挠度语义的 affine jet 处理，支持固定和浮动 基座。

递推参考示例覆盖固定基座多段模型、浮动基座跨段 tendon、浮动基座平面接触和 浮动基座 affine Ritz：

```shell
uv run --locked softarm build \
  examples/config/euler_bernoulli_pcs_recursive_n20.toml \
  --target matlab \
  --out examples/generated/euler_bernoulli_pcs_recursive_n20
uv run --locked softarm build \
  examples/config/extensible_euler_bernoulli_pcs_recursive_flying_n5.toml \
  --target matlab \
  --out examples/generated/extensible_euler_bernoulli_pcs_recursive_flying_n5
uv run --locked softarm build \
  examples/config/extensible_euler_bernoulli_pcs_recursive_flying_plane_contact_n1.toml \
  --target matlab \
  --out examples/generated/extensible_euler_bernoulli_pcs_recursive_flying_plane_contact_n1
uv run --locked softarm build \
  examples/config/euler_bernoulli_ritz_recursive_flying_n2.toml \
  --target matlab \
  --out examples/generated/euler_bernoulli_ritz_recursive_flying_n2
```

### 3.1 模型与积分

Extensible Euler–Bernoulli PCS 示例：

```toml
[model]
rod = "extensible_euler_bernoulli"
parameterization = "pcs"
segments = 2
inertia = "lumped"

[dynamics]
formulation = "symbolic_lagrange"

[integration]
method = "analytic"
```

`inertia` 可取 `distributed` 或 `lumped`。数值积分使用 `method = "gauss"` 并设置正整数 `order`。

Cosserat PCS 分布惯性要求显式选择至少二阶 Gauss 积分：

```toml
[model]
rod = "cosserat"
parameterization = "pcs"
segments = 1
inertia = "distributed"

[dynamics]
formulation = "symbolic_lagrange"

[integration]
method = "gauss"
order = 2
```

PAC 分布惯性要求至少四阶 Gauss–Legendre 积分；参考模型使用六阶：

```toml
[model]
rod = "euler_bernoulli"
parameterization = "pac"
segments = 2
inertia = "distributed"

[dynamics]
formulation = "symbolic_lagrange"

[integration]
method = "gauss"
order = 6
```

PAC 的仿射基固定，不使用额外配置表。集中惯性 PAC 可令 `inertia = "lumped"` 并省略 `[integration]`；分布惯性 PAC 不支持 `method = "analytic"`。

集中惯性使用 `inertia = "lumped"`，并省略 `[integration]`。每段六个坐标是 相对于可配置参考应变的增量，因此零坐标始终表示无应力参考构型。 当前参考包和完整数值验证覆盖单段 Cosserat PCS。配置和 builder 保留通用 `segments` 装配路径，但多段 Cosserat PCS 尚未验证，不作为当前版本的承诺能力。

Euler–Bernoulli Ritz 模型按

$$
\psi(\xi)=\sum_{k=0}^{m}c_k\xi^k
$$

配置两个弯曲平面的模态系数：

```toml
[model]
rod = "euler_bernoulli"
parameterization = "ritz"
segments = 2

[dynamics]
formulation = "symbolic_lagrange"

[ritz]
x = [0.0, 0.0, 1.5, -0.5]
y = [0.0, 0.0, 1.5, -0.5]
```

系数满足

$$
\psi(0)=\psi'(0)=0,\qquad \psi(1)=1.
$$

该配置可选择两种动力学装配。`recursive` 使用严格全局一阶 affine jet：基座 RPY 变换保持非线性，保留基座—软臂混合项，但不保留任意两个软臂坐标的乘积。

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

`kind = "unilateral"` 表示非负拉力单绳； `kind = "signed"` 表示拮抗绳索对的有符号等效通道。

`acceleration = "strict"` 生成绳长加速度约束求解器，适用于通道独立且 通道数不超过软臂坐标数的配置。`acceleration = "none"` 生成力映射， 适用于超驱动或相关通道组。

### 3.4 平面单点接触约束（实验性功能）

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

`plane_normal` 指向自由空间。该约束描述持续接触阶段；法向乘子的可行性输出 可用于上层接触状态管理。内置 `plane_point_contact` 同时支持两种动力学装配。 自定义约束通过 `register_constraint` 注册的 `ConstraintModel` 仅能供符号后端 使用；Recursive 后端目前只为内置平面点接触注册了工具点 jet 定义，其他约束族 会明确报告缺少 Recursive point-jet capability，不会回退到全局符号运动学。

## 4. 生成模型参考

### 4.1 模型清单

`manifest.json` 定义生成模型的公共元数据：

- `model`：Rod 理论、空间参数化、段数、基座模式、固定安装变换和中心线能力
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

`manifest.model.centerline` 为 `true` 的模型还提供 `softarm_centerline(q,p,xi)`，用于按归一化材料坐标采样各段中心线。

配置包含执行器时，生成包还会提供相应的 `softarm_actuator_*` 函数； 包含约束的 `symbolic_lagrange` 和 `recursive` 配置都会生成 `softarm_constraint_*` 函数。递推实现从工具点位置、Jacobian 和零加速度 bias 构造相同约束公式，不会为此展开全局符号质量矩阵或末端串联表达式。 Recursive bundle 还提供 `softarm_inverse_dynamics(q,dq,ddq,p)`；它不是 `symbolic_lagrange` MATLAB bundle 的公共函数。

### 4.3 CasADi bundle

安装固定版本的可选依赖，并从同一模型定义和局部 SymPy 公式生成 CasADi bundle：

```shell
uv sync --locked --extra casadi
uv run --locked softarm build examples/config/euler_bernoulli_pcs_n2.toml \
  --target casadi \
  --out build/euler_casadi
```

生成目录包含 `manifest.json`、`functions/*.casadi` 和普通数学手册 `softarm_model.tex`。载入接口不会依赖 acados：

```python
import casadi as ca
from softarm import load_casadi_bundle

bundle = load_casadi_bundle("build/euler_casadi")
f_expl = bundle.function("dynamics.generalized_force.explicit")

x = ca.MX.sym("x", f_expl.size1_in("x"))
u = ca.MX.sym("u", f_expl.size1_in("u"))
p = ca.MX.sym("p", f_expl.size1_in("p"))
expression = f_expl(x, u, p)
```

也可以直接使用 Python API：

```python
from softarm import generate_casadi_bundle

generate_casadi_bundle(
    plant,
    "build/model_casadi",
    actuation=actuation,
    constraint=constraint,
    affine_terms=None,
)
```

`constraint` 可为 `ConstraintConfig`，也可为已经物化的约束对象：Symbolic plant 接受 `ConstraintModel`，Recursive plant 接受 `ConstraintDefinition`。传入 `ConstraintConfig` 时生成器会自动选择对应表示。`p` 始终按 plant、actuator、 constraint 参数的顺序完整排列。

该表达式可由上层代码赋给 `AcadosModel.f_expl_expr`。每个适用的动力学变体 同时提供 `explicit` 和 `implicit` Function；具体 `x`、`xdot`、`u`、`z` 分块及 函数逻辑名称记录在 `manifest.dynamics`：

- `generalized_force`：软臂广义力、机体系基座扳手和世界系末端扳手。
- `tendon_force`：绳索张力、外部软臂广义力和两个外部扳手。
- `strict_tendon_acceleration`：严格绳坐标加速度；隐式形式以张力为代数量。
- `active_constraint`：已激活加速度约束；隐式形式以约束反力为代数量。

显式动力学使用 CasADi `solve` 求解质量矩阵或完整 KKT 系统。PCS、PAC 和 Cosserat 的递推隐式形式直接使用逆动力学残差，不构造质量矩阵；显式形式通过 对逆动力学关于广义加速度的 Jacobian 得到质量矩阵。Ritz 的全局一阶截断需要 专用 affine 变分装配，其质量矩阵和 bias 由 CasADi AD 从算法式 jet 构造。

PAC 中的仿射相位矩使用固定长度、可自动微分的级数递推。构建 PAC CasADi bundle 时必须显式指定 1 到 256 之间的项数，且该值会写入 manifest：

```shell
uv run --locked softarm build examples/config/euler_bernoulli_pac_distributed_n2.toml \
  --target casadi \
  --casadi-affine-terms 256 \
  --out build/pac_casadi
```

CasADi target 不接受 `--wolfram-cse`、`--wolfram-factor-terms` 或 `--tex-appendix`。单边张力和接触可行性应由上层优化约束表达，动力学 Function 本身不插入布尔断言。

Recursive 浮动 Ritz 与主动接触可直接生成 CasADi bundle：

```shell
uv run --locked softarm build examples/config/euler_bernoulli_ritz_recursive_flying_n2.toml \
  --target casadi --out build/ritz_recursive_casadi
uv run --locked softarm build \
  examples/config/extensible_euler_bernoulli_pcs_recursive_flying_plane_contact_n1.toml \
  --target casadi --out build/contact_recursive_casadi
```

后一个 bundle 可通过 `bundle.function("dynamics.active_constraint.explicit")`、`implicit` 和 `solution` 取得主动接触模型及 KKT 解函数。

### 4.4 LaTeX 数学手册

每次构建均生成 `softarm_model.tex`。`symbolic_lagrange` 文档依次给出模型摘要、 符号与参数、运动学、能量与动力学、执行器和约束，并与同一生成包中的数值函数 对应。`recursive` 文档记录递推算法概要；它不展开全局符号动力学表达式。

附加精确符号表达式：

```shell
uv run --locked softarm build examples/config/euler_bernoulli_ritz_n2.toml \
  --target matlab \
  --out build/euler \
  --tex-appendix
```

`--tex-appendix` 将本次构建所选规范化与 CSE 策略处理后的 $V,D,M,h,H_e,J_e,B_v$ 及适用的执行器、约束表达式加入附录。 该选项只适用于 `symbolic_lagrange`。

生成 PDF：

```shell
cd build/euler
pdflatex softarm_model.tex
```

TeX 文档依赖 `article`、`amsmath`、`amssymb`、`geometry` 和 `longtable`。 测试会先查找 `PATH` 中的 `pdflatex`；若未找到，则读取 `.softarm.local.toml` 的 `tools.pdflatex`。

## 5. MATLAB 接口

MATLAB 命令通过 `tools.matlab` 指定的可执行文件运行。例如重建 Simulink 模型：

```shell
uv run --locked softarm matlab -batch "addpath('matlab'); build_softarm_plant; build_softarm_constraint_models; build_softarm_actuator_models"
```

```matlab
addpath("matlab")
plant = softarm.loadModel("examples/generated/extensible_euler_bernoulli_pcs_lumped_n2");

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
    "examples/generated/extensible_euler_bernoulli_pcs_three_tendon_n2");
T = softarm.packActuatorInputs( ...
    plant.actuation,struct("t1",2,"t2",1.5,"t3",1));
[tau,isFeasible] = plant.actuation.force(q,T,plant.parameters);
ddq = plant.forwardDynamics( ...
    q,dq,tau,zeros(6,1),w,plant.parameters);
```

## 6. Simulink 模型

| 模型 | 用途 |
| --- | --- |
| `matlab/simulink/softarm_plant.slx` | 无约束 Plant；输入软臂广义力、基座扳手和末端扳手 |
| `matlab/simulink/softarm_constrained_plant.slx` | 约束 Plant；增加约束加速度并输出反力与可行性 |
| `matlab/simulink/softarm_actuator_force_block.slx` | 绳索拉力到软臂广义力的 Model Reference |
| `matlab/simulink/softarm_actuator_acceleration_block.slx` | 严格绳长加速度约束的 Model Reference |
| `examples/simulink/softarm_tendon_actuation_demo.slx` | 可切换绳索拉力/绳长加速度驱动示例 |
| `examples/simulink/softarm_flying_contact_demo.slx` | 浮动基座和平面单点接触示例 |

Plant 模型的 `Bundle` System Mask 用于选择生成包，并据此配置基座、软臂、 参数、执行器和约束维数。`backbone_poses` 端口按列堆叠每段末端的世界系 $4\times4$ 齐次变换，输出维数为 $16N\times1$。

在 Demo 的 `Plant` 块中，将 `Bundle directory` 设置为 MATLAB 字符串表达式， 即可切换模型，例如：

```matlab
'../generated/euler_bernoulli_ritz_two_signed_pairs_n2'
```

以上相对路径以 `examples/simulink` 中的 Demo 文件为基准。复用模型统一位于 `matlab/simulink`；Demo 的加载回调会自动将该目录加入 MATLAB 路径。

### 6.1 姿态回放

三个 Demo 将广义坐标记录到 `softarm_q_log`。仿真结束后运行：

```matlab
softarm_pose_playback
```

也可显式传入日志和生成包：

```matlab
softarm_pose_playback(softarm_q_log)
softarm_pose_playback(softarm_q_log, ...
    Bundle="examples/generated/euler_bernoulli_ritz_n2", ...
    SamplesPerSegment=16)
```

播放器根据 `q` 和生成包中的基座安装元数据重建软臂安装根与各段末端的世界系 位姿。`softarm_centerline` 在每段材料坐标上直接计算 PCS 或 Ritz 模型的中心线； `SamplesPerSegment` 默认为 16，可在 2 到 128 之间调节。外部模型未提供段内 几何时，播放器会警告并退回关键节点折线。

播放器用较长、较粗的 RGB 箭头显示刚体基座系 $H_{WB}$，并用短轴显示软臂安装 根 $H_{WB}H_{BM}$ 和各段末端姿态。世界原点是独立标记，不参与中心线连接。 输入必须是广义坐标日志，可采用 timeseries、timetable、structure-with-time、 包含 `softarm_q_log` 的 `Simulink.SimulationOutput`，或首列为时间的数值矩阵。 仅含段末端变换的 `backbone_poses` 日志不能用于离线回放。

## 7. 数理基础

### 7.1 符号与坐标约定

第 $i$ 段采用归一化材料坐标 $\xi\in[0,1]$，物理弧长坐标为 $s=L_i\xi$。各段广义坐标按段依次拼接为 $q$，速度为 $\dot q$。齐次变换

$$
{}^0H_i(q,\xi)=
\begin{bmatrix}
{}^0R_i(q,\xi)&{}^0r_i(q,\xi)\\
0&1
\end{bmatrix}
$$

把第 $i$ 段材料截面映射到 NED 基坐标系；多段模型由相邻段变换从基座到 末端依次左乘得到。

材料点的线速度和角速度 Jacobian 定义为

$$
J_{v,i}=\frac{\partial\,{}^0r_i}{\partial q},\qquad
J_{\omega,i}^{(:,j)}=
\mathrm{vex}\!\left(
\mathrm{skew}\!\left(
\frac{\partial\,{}^0R_i}{\partial q_j}\,{}^0R_i^T
\right)\right),
$$

其中 $\mathrm{skew}(A)=(A-A^T)/2$。线性 Euler 模型采用该式的一阶 小转角形式。末端空间 Jacobian 为

$$
J_e(q)=
\begin{bmatrix}
J_{v,e}J_{\omega,e}
\end{bmatrix}
\in\mathbb{R}^{6\times n_q}.
$$

### 7.2 PCS 运动学

每段 PCS 参数化使用常值角应变 $\kappa_i$、线应变 $\nu_i$ 和参考长度 $L_i$。 令

$$
\Omega_i(\xi)=L_i\xi\widehat{\kappa_i},\qquad
z_i=(L_i\xi)^2\kappa_i^T\kappa_i,
$$

则单段位姿由统一的 SE(3) 指数映射给出：

$$
R_i=I+S(z_i)\Omega_i+C(z_i)\Omega_i^2,
$$

$$
r_i=\left[I+C(z_i)\Omega_i+T(z_i)\Omega_i^2\right]L_i\xi\nu_i,
$$

其中 $S$、$C$ 和 $T$ 均在零曲率处采用解析延拓。Euler–Bernoulli PCS 使用

$$
q_i=[b_{x,i},b_{y,i}]^T,\quad
\kappa_i=[-b_{y,i}/L_i,b_{x,i}/L_i,0]^T,\quad
\nu_i=[0,0,1]^T,
$$

因此广义坐标中没有长度项。Extensible Euler–Bernoulli PCS 使用

$$
q_i=[b_{x,i},b_{y,i},l_i]^T,\quad
\kappa_i=[-b_{y,i}/L_{0,i},b_{x,i}/L_{0,i},0]^T,\quad
\nu_i=[0,0,l_i/L_{0,i}]^T.
$$

### 7.3 PAC 运动学

PAC 每段采用三维仿射曲率坐标。Euler–Bernoulli 段为 $q_i=[c_{0,i},c_{1,i},\phi_i]^T$，Extensible Euler–Bernoulli 段再加入绝对长度 $l_i$。令

$$
\alpha_i(\xi)=c_{0,i}\xi+\frac12c_{1,i}\xi^2,
$$

$$
\mathcal C_n=\int_0^\xi v^n\cos\left(c_{0,i}v+\frac12c_{1,i}v^2\right)dv,
\qquad
\mathcal S_n=\int_0^\xi v^n\sin\left(c_{0,i}v+\frac12c_{1,i}v^2\right)dv.
$$

局部变换为

$$
R_i(\xi)=R_z(\phi_i)R_y(\alpha_i(\xi)),
\qquad
r_i(\xi)=\ell_i
\begin{bmatrix}
\cos\phi_i\,\mathcal S_0\\
\sin\phi_i\,\mathcal S_0\\
\mathcal C_0
\end{bmatrix},
$$

其中 Euler–Bernoulli 取 $\ell_i=L_i$，Extensible Euler–Bernoulli 取 $\ell_i=l_i$。因此前者严格满足 $\|\partial r_i/\partial\xi\|=L_i$，后者的中心线弧长为 $l_i$。 积分特殊函数在 $c_1=0$ 和完全直杆处使用解析延拓。

以参考弧长为材料坐标，段内连续应变为

$$
\begin{aligned}
\text{Euler–Bernoulli:}\quad
&\kappa_i(\xi)=
\begin{bmatrix}0&(c_{0,i}+c_{1,i}\xi)/L_i&0\end{bmatrix}^T,
&&\nu_i=\boldsymbol e_3,\\
\text{Extensible Euler--Bernoulli:}\quad
&\kappa_i(\xi)=
\begin{bmatrix}0&(c_{0,i}+c_{1,i}\xi)/L_{0,i}&0\end{bmatrix}^T,
&&\nu_i=(l_i/L_{0,i})\boldsymbol e_3,
\end{aligned}
$$

其中 $\boldsymbol e_3=[0,0,1]^T$。段起点截面旋转满足 $R_i(0)=R_z(\phi_i)$；$\phi_i$ 表示相邻段之间的离散材料截面转角，不对应段内的 连续扭转应变。直杆中心线对 $\phi_i$ 不敏感，多方向绳索 Jacobian 可能在该构型降秩。

令 $c_i=[c_{0,i},c_{1,i}]^T$、 $H=\begin{bmatrix}1&1/2\\1/2&1/3\end{bmatrix}$。各向异性 Euler PAC 弹性能为

$$
V_{e,i}=\frac12\frac{EI_{y,i}\cos^2\phi_i+EI_{x,i}\sin^2\phi_i}{L_i}
c_i^THc_i+\frac12\frac{GJ_i}{L_i}\phi_i^2.
$$

Extensible Euler–Bernoulli PAC 弹性能为

$$
V_{e,i}=\frac12(k_{bx,i}\cos^2\phi_i+k_{by,i}\sin^2\phi_i)c_i^THc_i
+\frac12k_{\phi,i}\phi_i^2+\frac12k_{l,i}(l_i-L_{0,i})^2.
$$

### 7.4 Euler–Bernoulli 与 Extensible Euler–Bernoulli Ritz 运动学

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

该模型适用于小挠度、小斜率和小转角工况。Euler–Bernoulli Ritz 令线性轴向应变为零， 但其中心线导数满足

$$
\left\|\frac{\partial r}{\partial\xi}\right\|=
\sqrt{L^2+a_x^2(\psi_x')^2+a_y^2(\psi_y')^2},
$$

因此不精确保持弧长；相对 $L$ 的几何伸长为斜率的二阶量，被一阶模型忽略。 Extensible Euler–Bernoulli Ritz 再加入 $w_z(\xi)=a_z\psi_z(\xi)$，并保留线性轴向应变 $\varepsilon_z=(a_z/L)\psi_z'(\xi)$；横向斜率引起的二阶轴向应变同样不计。

Euler 段的弯曲势能为

$$
V_{e,i}=
\frac12\frac{EI_{y,i}a_{x,i}^2}{L_i^3}
\int_0^1\left(\psi_x''\right)^2d\xi
+
\frac12\frac{EI_{x,i}a_{y,i}^2}{L_i^3}
\int_0^1\left(\psi_y''\right)^2d\xi.
$$

### 7.5 Cosserat PCS 运动学

每段使用参考应变 $\kappa_{0,i},\nu_{0,i}$ 和六维应变增量：

$$
\kappa_i=\kappa_{0,i}+\delta\kappa_i,\qquad
\nu_i=\nu_{0,i}+\delta\nu_i.
$$

令 $\Omega_i(\xi)=L_i\xi\widehat{\kappa_i}$、 $z_i(\xi)=(L_i\xi)^2\kappa_i^T\kappa_i$，并定义

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
H_i(\xi)=\begin{bmatrix}R_i&r_i\\
0&1\end{bmatrix}.
$$

参考应变默认为 $\kappa_0=0,\nu_0=[0,0,1]^T$，也可通过运行时参数设置 预弯、预扭或参考剪切。弹性能为

$$
V_{e,i}=\frac{L_i}{2}\delta\xi_i^T
\mathrm{diag}(EI_x,EI_y,GJ,GA_x,GA_y,EA)\delta\xi_i.
$$

### 7.6 能量、质量矩阵与偏置力

对于均匀分布惯性，$m_i$ 和 $I_i$ 表示整段总质量和整段局部惯量。 归一化材料坐标下 $dm=m_i\,d\xi$，每段质量矩阵贡献为

$$
M_i(q)=\int_0^1\left[
m_iJ_{v,i}^TJ_{v,i}
+J_{\omega,i}^T\,{}^0R_iI_i{}^0R_i^T\,J_{\omega,i}
\right]d\xi.
$$

线性 Ritz 模型采用参考姿态 ${}^0R_i(0)$ 的惯量旋转，以保持运动学一阶、 能量二阶的一致截断。PCS `lumped` 模式在 $\xi=1/2$ 处计算该段贡献。 末端负载贡献为

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

PCS `lumped` 模式的对应项为 $-m_i g z_i(q,1/2)$。Extensible Euler–Bernoulli PCS 的有效弹性势能为

$$
V_{e,i}^{\mathrm{EEB,PCS}}=
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

`symbolic_lagrange` 直接构造上述全局符号 $M$ 和 $h$。PCS、PAC 和 Cosserat PCS 的 Exact-SE(3) `recursive` 路径定义逆动力学映射 $r(q,\dot q,\ddot q)$：先逐段前向传播速度与加速度，再反向回拉空间力旋量。 MATLAB 代码按

$$
h(q,\dot q)=r(q,\dot q,0),\qquad
M_{:j}(q)=r(q,0,e_j)-r(q,0,0)
$$

数值装配质量矩阵和偏置项，其中 $e_j$ 是第 $j$ 个单位加速度。CasADi 隐式形式 直接使用 $r-Q$，显式形式通过 $\partial r/\partial\ddot q$ 得到质量矩阵。

displacement-Ritz 不使用上述空间力递推。它传播全局一阶 affine 变换 jet，以固定 Gauss 求积直接累加分布质量、惯量、重力、弹性和阻尼；浮动基座变换保持完整 非线性。MATLAB 由该 jet 形成 $M,h$，CasADi 通过算法式组合和 AD 形成相同量。 所有显式自由动力学和 KKT 动力学都求解线性系统，不在生成代码中形成显式逆矩阵。 两种装配方法保持相同的 MATLAB 公共动力学接口，Recursive 路径不构造随总段数 膨胀的全局 SymPy $M,h$。

设 $S_a$ 为软臂坐标选择矩阵，$J_B$ 为基座刚体 Jacobian， $R_{WB}$ 为机体系到 NED 世界系的旋转，则

$$
Q=S_a\tau_a+B_v(q)w_B+J_e^Tw_e,\qquad
B_v=J_B^T\mathrm{diag}(R_{WB},R_{WB}).
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

`softarm.linearize` 在指定工作点对该状态方程进行中心差分，返回 $A=\partial\dot x/\partial x$ 以及软臂广义力、基座扳手和末端扳手对应的 输入矩阵。

### 7.7 材料坐标积分

解析积分为默认策略。Gauss–Legendre 积分显式配置为 $n$ 阶时，

$$
\int_0^1f(\xi)d\xi\approx
\frac12\sum_{k=1}^{n}w_k
f\!\left(\frac{\eta_k+1}{2}\right),
$$

其中 $(\eta_k,w_k)$ 为区间 $[-1,1]$ 上的 Gauss–Legendre 节点和权重。 构建过程严格采用配置指定的积分策略。

### 7.8 绳索驱动

第 $a$ 个 Extensible Euler–Bernoulli PCS 单绳通道的长度坐标为

$$
y_a=\sum_{i\in\mathcal P_a}\left[
l_i-r_{ai}
\left(
b_{x,i}\cos\alpha_{ai}+b_{y,i}\sin\alpha_{ai}
\right)
\right].
$$

`signed` 通道采用相对截面中心对称布置的拮抗绳索半差动长度， 其轴向分量为零。Euler 通道使用 Ritz 端部斜率构造相应弯曲项。

PAC 段使用总转角 $\bar c_i=c_{0,i}+c_{1,i}/2$。截面极角为 $\theta_{ai}$ 的单绳 span 长度为

$$
y_{ai}=\ell_i-r_{ai}\bar c_i\cos(\theta_{ai}-\phi_i),
$$

其中 Euler 取 $\ell_i=L_i$，Extensible Euler–Bernoulli 取 $\ell_i=l_i$； `signed` 通道去掉 $\ell_i$。严格加速度模式仍要求无应力参考构型下 $J_a$ 满行秩，降秩配置直接拒绝且不做正则化。

Cosserat PCS 在偏置 $r=[r\cos\theta,r\sin\theta,0]^T$ 处使用精确常应变绳长

$$
\ell_+=L\|\nu+\kappa\times r\|,
$$

单绳通道取 $\ell_+$，signed 对取 $\frac{L}{2}(\|\nu+\kappa\times r\|-\|\nu-\kappa\times r\|)$。

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
\ddot q\\
T
\end{bmatrix}
=\begin{bmatrix}
Q-h\\
\ddot y_{\mathrm{cmd}}-\dot J_a\dot q_a
\end{bmatrix}.
$$

### 7.9 加速度级约束

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
\ddot q\\
\lambda
\end{bmatrix}
=\begin{bmatrix}
Q-h\\
b-\dot A\dot q
\end{bmatrix}.
$$

`plane_point_contact` 使用单边法向约束，并以平滑库仑模型将切向摩擦合入 反力映射。求解结果同时提供法向反力可行性和 KKT 条件诊断。

## 8. 示例索引

表中标识符同时对应 `examples/config/<标识符>.toml` 和 `examples/generated/<标识符>/`。

| 示例 | 标识符 | 动力学装配 |
| --- | --- | --- |
| 两段 Extensible Euler–Bernoulli PCS，中点集中惯性 | `extensible_euler_bernoulli_pcs_lumped_n2` | `symbolic_lagrange` |
| 两段 Extensible Euler–Bernoulli PCS，分布惯性 | `extensible_euler_bernoulli_pcs_distributed_n2` | `symbolic_lagrange` |
| 两段 Euler–Bernoulli PAC，分布惯性 | `euler_bernoulli_pac_distributed_n2` | `symbolic_lagrange` |
| 两段 Extensible Euler–Bernoulli PAC，分布惯性 | `extensible_euler_bernoulli_pac_distributed_n2` | `symbolic_lagrange` |
| 两段 Euler–Bernoulli Ritz | `euler_bernoulli_ritz_n2` | `symbolic_lagrange` |
| 两段 Euler–Bernoulli PCS | `euler_bernoulli_pcs_n2` | `symbolic_lagrange` |
| 两段 Extensible Euler–Bernoulli Ritz | `extensible_euler_bernoulli_ritz_n2` | `symbolic_lagrange` |
| 两段 Extensible Euler–Bernoulli PCS，三绳驱动 | `extensible_euler_bernoulli_pcs_three_tendon_n2` | `symbolic_lagrange` |
| 两段 Extensible Euler–Bernoulli PCS，signed 绳索对 | `extensible_euler_bernoulli_pcs_signed_pair_n2` | `symbolic_lagrange` |
| 两段 Euler–Bernoulli Ritz，两个 signed 绳索对 | `euler_bernoulli_ritz_two_signed_pairs_n2` | `symbolic_lagrange` |
| 浮动基座 Extensible Euler–Bernoulli PCS，平面单点接触 | `extensible_euler_bernoulli_pcs_flying_plane_contact_n1` | `symbolic_lagrange` |
| 浮动基座 Extensible Euler–Bernoulli PCS，递推平面单点接触 | `extensible_euler_bernoulli_pcs_recursive_flying_plane_contact_n1` | `recursive` |
| 浮动基座两段 Euler–Bernoulli Ritz | `euler_bernoulli_ritz_flying_n2` | `symbolic_lagrange` |
| 浮动基座两段 Euler–Bernoulli affine Ritz | `euler_bernoulli_ritz_recursive_flying_n2` | `recursive` |
| 单段 Cosserat PCS，分布惯性与三绳驱动 | `cosserat_pcs_distributed_tendon_n1` | `symbolic_lagrange` |
| 单段 Cosserat PCS，中点集中惯性 | `cosserat_pcs_lumped_n1` | `symbolic_lagrange` |
| 二十段 Euler–Bernoulli PCS，固定基座、第一段三绳驱动 | `euler_bernoulli_pcs_recursive_n20` | `recursive` |
| 五段 Extensible Euler–Bernoulli PCS，浮动基座、跨段 signed 绳索对 | `extensible_euler_bernoulli_pcs_recursive_flying_n5` | `recursive` |

MATLAB 绳索驱动脚本见 [`examples/actuation`](examples/actuation/)。
