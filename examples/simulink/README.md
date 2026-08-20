# Simulink Demo

本目录包含可直接运行的 SoftArm 顶层示例：

- `softarm_tendon_actuation_demo.slx`：可通过 `Actuator mode` Mask 切换三绳拉力或严格绳长加速度驱动。
- `softarm_flying_contact_demo.slx`：浮动基座平面接触。

Demo 引用的通用 Plant 和执行器模型位于 `../../matlab/simulink/`，模型加载回调
会自动加入所需 MATLAB 路径。默认生成包位于 `../generated/`。

所有 SLX 均由 `matlab/build_softarm_*.m` 重建，不应手工维护二进制模型。
