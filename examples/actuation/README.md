# 通用绳索驱动示例

本目录演示 MATLAB 中的 `softarm_actuator_*` 生成接口。执行器模块将绳索输入
映射为 Plant 的软臂广义力 `tau`。

- `run_three_tendon_force.m`：三根相隔 120° 的非负拉力驱动两段 Extensible Kirchhoff PCS rod。
- `run_three_tendon_acceleration.m`：严格约束三根贯穿绳索的长度加速度，并输出求解拉力与可行性。
- `run_signed_force.m`：两个理想成对通道接受正负等效拉力，其 Jacobian 的轴向列为零。

重新生成参考包：

```powershell
softarm build examples/config/extensible_kirchhoff_pcs_three_tendon_n2.toml --out examples/generated/extensible_kirchhoff_pcs_three_tendon_n2
softarm build examples/config/extensible_kirchhoff_pcs_signed_pair_n2.toml --out examples/generated/extensible_kirchhoff_pcs_signed_pair_n2
softarm build examples/config/euler_bernoulli_ritz_two_signed_pairs_n2.toml --out examples/generated/euler_bernoulli_ritz_two_signed_pairs_n2
```

正拉力遵循 `tau = -Ja.'*tension`。`signed` 的负值表示拮抗绳索对的反向等效
作用；各物理绳索仍以非负拉力工作。

## Simulink 示例

打开 `examples/simulink/softarm_tendon_force_demo.slx` 或
`examples/simulink/softarm_tendon_acceleration_demo.slx`。模型启动回调将加载
`matlab` 与 `matlab/simulink` 路径和默认的
`extensible_kirchhoff_pcs_three_tendon_n2` 生成包。

拉力 Demo 中双击 `Tension command` 修改三维拉力向量；加速度 Demo 中双击
`Tendon acceleration` 修改三维绳长加速度。点击工具栏 Run，通过
`State scope` 查看广义坐标与速度。

Bundle 选择位于 `softarm_plant.slx` 定义的 System Mask 中。双击 Demo 的
`Plant` 块，将 `Bundle directory` 设置为以下模型参数表达式并点击 Apply，
即可切换到两组 signed 绳索对驱动的 Euler 模型：

```matlab
'../generated/euler_bernoulli_ritz_two_signed_pairs_n2'
```

该字段采用 MATLAB 字符串表达式，因此路径使用单引号。切换后输入命令调整为
二维，Euler 的 `q/dq` 各为四维，System Mask 同步更新 Plant 与执行器维数。

在自定义顶层模型中复用时，添加 Model block，选择
`softarm_actuator_force_block` 或 `softarm_actuator_acceleration_block`，
并按示例连接到 `softarm_plant`。这些可复用模型位于 `matlab/simulink`。
