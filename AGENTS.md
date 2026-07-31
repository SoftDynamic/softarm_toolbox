# SoftArm Toolbox 维护指南

本文件适用于仓库根目录及全部子目录，记录模型实现、代码生成和验证过程中的
内部约束。面向使用者的模型配置、接口和数学定义维护在 `README.md`。

## 1. 架构边界

- SymPy 表达式是运动学、能量、动力学、执行器和约束公式的规范表示，
  符号推导与材料坐标积分均在该层完成。
- Wolfram 后端提供可选的表达式优化与 CSE；MATLAB 和 Simulink 消费生成的
  数值函数。
- 同一公式不得在 Python、Wolfram 和 MATLAB 层分别维护。
- `SymbolicPlant` 的公共结果保持为 `mass`、`potential`、`damping`、
  `kinematics` 和 `end_jacobian`；`bias` 由统一能量公式构造。
- 结构参数进入 TOML。仿真期间需要调节的物理参数进入有序运行时参数向量 `p`。
- 本机工具路径写入 `.softarm.local.toml`，该文件属于机器本地配置。

## 2. 模型实现

### 2.1 修改现有模型

1. 单段运动学在 `src/softarm/geometry.py` 中修改，公式使用 SymPy 对象构造。
2. 参数、能量和惯性装配在 `src/softarm/derive.py` 中修改。
3. 配置字段及校验在 `src/softarm/config.py` 中修改，并同步更新
   `examples/config/`。
4. 生成函数签名或输出维度发生变化时，同步更新
   `src/softarm/codegen/matlab.py`、MATLAB 加载层和
   `matlab/build_softarm_plant.m`。
5. 重新生成受影响的 `examples/generated/` 参考包，并执行完整验证。

### 2.2 增加几何或梁模型

- Builder 接收 `ModelConfig` 并返回 `SymbolicPlant`。
- 仓库内置模型注册到 `derive.py` 的 `_MODEL_BUILDERS`。
- 外部 Python 模型在调用 `derive` 前通过
  `softarm.register_model(name, builder)` 注册。CLI 集成需要在启动阶段显式导入
  注册模块。
- 新模型应明确选择适用假设，包括小挠度或大转角、可伸长或不可伸长、
  有剪切或无剪切。不同阶次或相互冲突的假设应定义为独立模型。
- 通用能量装配应复用 `derive.py` 中的公共路径；结构常量应通过
  `ModelConfig` 和 TOML 显式表达。

### 2.3 特殊函数与可去奇点

- `SincSqrt` 和 `CoscSqrt` 表示 PCC 零曲率处的解析延拓，并提供一、二阶解析导数。
- MATLAB 数值函数在 $|\rho^2|<10^{-8}$ 区间使用 Taylor 多项式。
- 新增特殊函数时，在 `src/softarm/special.py` 中定义 SymPy 函数、解析导数和
  独立数值实现。
- 同一 AST 节点的往返和输出规则应同步加入 `src/softarm/ast.py`、
  `src/softarm/backends/wolfram_bridge.wls` 和 MATLAB printer。
- 奇点处理应保持状态变量和动力学矩阵的原始定义，禁止通过修改状态的 epsilon
  或人为质量正则项替代解析延拓。

## 3. 积分与后端

- `integration.method = "analytic"` 表示由 SymPy 执行材料坐标解析积分。
- `integration.method = "gauss"` 表示显式生成指定阶数的 Gauss–Legendre 求和。
- 积分方法由配置确定；解析积分错误应作为构建错误报告。
- SymPy 和 Wolfram 后端必须对同一中间表达式保持数值一致。

## 4. 执行器与约束

### 4.1 执行器

- Plant 公共输入为软臂广义力、机体系基座扳手和世界系末端扳手。
- 绳索、气腔或电机模型实现 `ActuatorMap`，并在 Plant 外部计算
  $\tau_a=B(q_a,p_u)u$。
- 自定义执行器通过 `softarm.register_actuator(name,builder)` 注册；
  builder 返回完全由 SymPy 表达式构成的 `ActuationModel`。
- 严格绳长加速度模式使用无正则化 KKT 系统，仅接受线性独立且数量不超过
  软臂坐标数的通道。

### 4.2 约束

- `ConstraintModel` 描述约束值、Jacobian、速度偏置、反力映射和稳定化参数。
- 理想约束使用 $G=A^T$；非理想反力可提供独立的 $G$。
- 内置约束求解器使用无正则化 KKT 方程。
- `plane_point_contact` 表示已激活的持续接触阶段。接近、碰撞和脱离由上层
  模式管理器根据间隙和反力可行性管理。

## 5. 生成包与数学手册

- `manifest.json` 是生成包的元数据来源，包含模型、分组坐标、有序参数、
  执行器和约束描述。
- MATLAB 加载层以当前 manifest 结构为准。修改结构时，应同步更新加载器、
  参考包和测试。
- 公共 MATLAB 函数签名保持稳定；签名调整属于跨 Python 生成器、MATLAB
  加载层和 Simulink 模型的联合变更。
- `softarm_model.tex` 与 MATLAB 代码共享同一 `SymbolicPlant`、
  `ActuationModel` 和 `ConstraintModel`。
- 数学手册正文使用命名装配公式。`--tex-appendix` 用于输出优化和 CSE 后的
  精确表达式；大型浮动基座模型的附录体量应在提交前评估。

## 6. MATLAB 与 Simulink 维护

- SLX 文件由 MATLAB 构建脚本生成，禁止将手工编辑作为模型维护路径。
- 重建顺序为：

  ```matlab
  build_softarm_plant
  build_softarm_constraint_models
  build_softarm_actuator_models
  ```

- `softarm_plant.slx` 和 `softarm_constrained_plant.slx` 的已有端口编号保持稳定。
- `backbone_poses` 按列堆叠每段末端的 $4\times4$ 世界系齐次变换，
  维数为 $16N\times1$。
- Demo 顶层使用 `To Workspace` 记录 `q`，变量名为 `softarm_q_log`，
  `Sample time` 设为 `-1`。
- 姿态重建与绘图在仿真结束后由 `softarm_pose_playback` 执行。
- Model Reference 反馈路径中的 Memory 块用于打断 Simulink 的保守代数环判断。
- Demo 的 Simulation Pacing 保持关闭。

## 7. 验证

Python 测试：

```shell
python -m pytest
```

MATLAB 与 Simulink 测试：

```shell
matlab -batch "addpath('matlab'); r=runtests('matlab/tests'); assertSuccess(r)"
```

模型变更至少覆盖：

- 零构型或参考构型
- 运动学有限差分
- 质量矩阵对称性与名义正定性
- 独立能量或积分计算
- Coriolis 能量恒等式
- 外力虚功
- PCC 零曲率处的函数值及一、二阶导数
- SymPy、Wolfram 和 MATLAB 数值一致性
- 执行器虚功、轴向模态和严格加速度可行性
- 约束 Jacobian、反力映射和摩擦耗散
- 单段与多段组合
- manifest、MATLAB 加载、数值线性化和 Simulink 短时运行
