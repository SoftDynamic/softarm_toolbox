# 通用绳索驱动示例

这里的 MATLAB 脚本使用同一组 `softarm_actuator_*` 生成接口。Plant 本身仍只接收广义力 `tau`；绳索映射作为 Plant 外部模块使用。

- `run_three_tendon_force.m`：三根相隔 120° 的非负拉力驱动两段可伸缩 PCC。
- `run_three_tendon_acceleration.m`：严格约束三根贯穿绳索的长度加速度，并输出求解拉力与可行性。
- `run_signed_force.m`：两个理想成对通道接受正负等效拉力，其 Jacobian 的 PCC 轴向列为零。

重新生成参考包：

```powershell
softarm build examples/config/pcc_three_tendon_extensible_n2.toml --out examples/generated/pcc_three_tendon_extensible_n2
softarm build examples/config/pcc_signed_pair_n2.toml --out examples/generated/pcc_signed_pair_n2
softarm build examples/config/euler_two_signed_pairs_n2.toml --out examples/generated/euler_two_signed_pairs_n2
```

正拉力遵循 `tau = -Ja.'*tension`。`signed` 的负值表示等效配对的反向作用，不是单根物理绳承受压力。

## 无代码 Simulink 运行

直接打开仓库根目录的 `softarm_tendon_force_demo.slx` 或 `softarm_tendon_acceleration_demo.slx`。模型会自动选择 `pcc_three_tendon_extensible_n2`，无需执行 `addpath` 或初始化脚本。

拉力 Demo 中双击 `Tension command` 修改三维拉力向量；加速度 Demo 中双击 `Tendon acceleration` 修改三维绳长加速度。点击工具栏 Run，通过 `State scope` 查看广义坐标与速度。

Bundle 选择位于 `softarm_plant.slx` 定义的 System Mask 中。双击 Demo 里的 `Plant` 块，把 `Bundle directory` 从默认值改为下面的模型参数表达式并点击 Apply，即可切换到“两对 tendon + Euler”：

```matlab
'examples/generated/euler_two_signed_pairs_n2'
```

路径外的单引号是必需的。切换后输入命令会自动成为二维，Euler 的 `q/dq` 各为四维；Plant 和驱动器之间的连线无需修改。

在自己的顶层模型中复用时，添加 Model block 并分别选择 `softarm_actuator_force_block` 或 `softarm_actuator_acceleration_block`，再按示例图连接到 `softarm_plant`。
