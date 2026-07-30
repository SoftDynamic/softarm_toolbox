clear; close all;
addpath("functions");

%% 1. 参数定义
N = 2;                                  % 软臂节数
g = 9.81;                               % NED 坐标系重力加速度 (m/s^2)

M_i = 0.2;                              % 每节软臂质量 (kg)
L_i_0 = 0.5;                            % 每节软臂固定长度 (m)
I_i_local = diag([0.002, 0.002, 1e-4]); % 每节软臂的局部惯性矩 (kg*m^2)
R = 0.0275;                             % 驱动绳分布半径 (m)
Me = 0.1;                               % 末端负载质量 (kg)
Ie_local = diag([1e-2, 1e-2, 1e-2]);    % 末端负载的局部惯性矩 (kg*m^2)

EI = 1.2;                               % Euler-Bernoulli 弯曲刚度 (N*m^2)
d_x = 0.05;                             % x 方向弯曲阻尼
d_y = 0.05;                             % y 方向弯曲阻尼

%% 2. 符号变量与 displacement 状态
x = sym('x', [1 N], 'real');
y = sym('y', [1 N], 'real');
dotX = sym('dotX', [1 N], 'real');
dotY = sym('dotY', [1 N], 'real');

u_acc = sym('u_acc', [2 1], 'real');

% 末端外力仅保留三维力，不引入外力矩接口
syms fx fy fz real
F_ext = [fx; fy; fz];

% q = [x1; y1; x2; y2]，状态 X = [q; dotQ]
q = reshape([x; y], [], 1);
dotQ = reshape([dotX; dotY], [], 1);

%% 3. 运动学分析
H_0 = sym(zeros(4, 4, N));
H_0_com = sym(zeros(4, 4, N));
Jv_com = sym(zeros(3, 2*N, N));
Jw = sym(zeros(3, 2*N, N));
Jw_com = sym(zeros(3, 2*N, N));

for i = 1:N
    H_i = euler_homogeneous(x(i), y(i), L_i_0);
    H_com_i_local = euler_homogeneous(x(i)/2, y(i)/2, L_i_0/2);

    if i == 1
        H_0(:,:,i) = H_i;
        H_0_com(:,:,i) = H_com_i_local;
        R_previous = sym(eye(3));
        Jw_previous = sym(zeros(3,2*N));
    else
        H_0(:,:,i) = H_0(:,:,i-1) * H_i;
        H_0_com(:,:,i) = H_0(:,:,i-1) * H_com_i_local;
        R_previous = H_0(1:3,1:3,i-1);
        Jw_previous = Jw(:,:,i-1);
    end

    % 直接由旋转矩阵导数构造当前节的局部角速度雅可比
    R_i = H_i(1:3,1:3);
    Jw_i_local = [vex3(diff(R_i, x(i)) * R_i'), ...
                  vex3(diff(R_i, y(i)) * R_i')];

    % 质心旋转使用半弯曲量，与质心齐次变换一致
    R_com_i_local = H_com_i_local(1:3,1:3);
    Jw_com_i_local = [vex3(diff(R_com_i_local, x(i)) * R_com_i_local'), ...
                      vex3(diff(R_com_i_local, y(i)) * R_com_i_local')];

    current_columns = (2*i-1):(2*i);
    Jw(:,:,i) = Jw_previous;
    Jw(:,current_columns,i) = R_previous * Jw_i_local;
    Jw_com(:,:,i) = Jw_previous;
    Jw_com(:,current_columns,i) = R_previous * Jw_com_i_local;
    Jv_com(:,:,i) = jacobian(H_0_com(1:3,4,i), q);
end

p_end = H_0(1:3,4,N);
Jv_end = jacobian(p_end, q);
Jw_end = Jw(:,:,N);

%% 4. 绳索 displacement 映射
Lx_total = 0;
Ly_total = 0;
for i = 1:N
    Lx_total = Lx_total + R * x(i);
    Ly_total = Ly_total + R * y(i);
end
Ja = jacobian([Lx_total; Ly_total], q);
dotJa_dotQ = jacobian(Ja * dotQ, q) * dotQ;

%% 5. 能量模型
% 5.1 软臂各节的平动与转动动能
K_arm = 0;
M_arm = sym(zeros(2*N));
for i = 1:N
    M_trans_i = M_i * (Jv_com(:,:,i)' * Jv_com(:,:,i));
    K_trans = 0.5 * (dotQ' * M_trans_i * dotQ);

    R_com_i = H_0_com(1:3,1:3,i);
    I_i_base = R_com_i * I_i_local * R_com_i';
    M_rot_i = Jw_com(:,:,i)' * I_i_base * Jw_com(:,:,i);
    K_rot = 0.5 * (dotQ' * M_rot_i * dotQ);

    K_arm = K_arm + K_trans + K_rot;
    M_arm = M_arm + M_trans_i + M_rot_i;
end

% 5.2 末端负载的平动与转动动能
R_end = H_0(1:3,1:3,N);
Ie_base = R_end * Ie_local * R_end';
M_load = Me * (Jv_end' * Jv_end) + Jw_end' * Ie_base * Jw_end;
K_load = 0.5 * (dotQ' * M_load * dotQ);
K_sum = K_arm + K_load;

% 5.3 NED 重力势能与 Euler-Bernoulli 弯曲势能
P_gravity = 0;
for i = 1:N
    P_gravity = P_gravity - M_i * g * H_0_com(3,4,i);
end
P_gravity = P_gravity - Me * g * p_end(3);

k_bending = EI / L_i_0;
P_elastic = 0;
for i = 1:N
    P_elastic = P_elastic + ...
        0.5 * k_bending * (x(i)^2 + y(i)^2);
end
P_sum = P_gravity + P_elastic;

%% 6. 动力学矩阵提取
disp('正在计算质量矩阵...');
% K_sum = 0.5*dotQ'*Mass_Matrix*dotQ，直接累加各刚体项可避免
% 对展开后的标量动能进行两次代价高昂的符号求导。
Mass_Matrix = M_arm + M_load;

disp('正在计算保守力 (G+K)...');
G_K_Matrix = jacobian(P_sum, q)';

D_Matrix = diag(repmat([d_x, d_y], 1, N)) * dotQ;

disp('正在计算科里奥利项...');
% Christoffel 向量的等价矢量形式，避免构造并化简完整 C 矩阵。
C_vector = jacobian(Mass_Matrix * dotQ, q) * dotQ - ...
           0.5 * jacobian(dotQ' * Mass_Matrix * dotQ, q)';

% 总偏置力
n_vec = G_K_Matrix + D_Matrix + C_vector;
%% 7. 导出仿真函数
disp('正在导出函数...');
if ~isfolder('output')
    mkdir('output');
end

matlabFunction(Mass_Matrix, ...
    'File', 'output/get_MassMatrix_Euler', ...
    'Vars', {q}, ...
    'Optimize', true);

matlabFunction(n_vec, ...
    'File', 'output/get_BiasForce_Euler', ...
    'Vars', {q, dotQ}, ...
    'Optimize', false);

matlabFunction(Jv_end, ...
    'File', 'output/get_EndJacobian_Euler', ...
    'Vars', {q}, ...
    'Optimize', false);

matlabFunction(H_0, ...
    'File', 'output/get_homogeneous_Euler', ...
    'Vars', {q}, ...
    'Optimize', false);

%% 8. displacement 约束与状态空间模型
% M(q)ddq + n(q,dq) = Ja'F_rope + Jv_end'F_ext
% Ja*ddq + dot(Ja)*dotQ = u_acc

if N <= 2
    
    LHS = [Mass_Matrix, -Ja';
           Ja,           zeros(2,2)];
    RHS = [Jv_end' * F_ext - n_vec;
           u_acc - dotJa_dotQ];
    
    % 对 KKT 增广方程作等价的零空间消元。B 是 Ja 的右逆，Z 张成
    % null(Ja)，因此绳索约束力在 Z' 投影后消失，只需求解 2x2 系统。
    B_constraint = Ja' / (Ja * Ja');
    Z_constraint = sym(zeros(2*N, 2*(N-1)));
    for i = 1:N-1
        Z_constraint(2*i-1,2*i-1) = 1;
        Z_constraint(2*i,2*i) = 1;
        Z_constraint(2*N-1,2*i-1) = -1;
        Z_constraint(2*N,2*i) = -1;
    end
    
    if N ==2
        M_reduced = Z_constraint' * Mass_Matrix * Z_constraint;
        bias_reduced_rhs = -Z_constraint' * n_vec + ...
            Z_constraint' * Mass_Matrix * B_constraint * dotJa_dotQ;
        input_reduced_rhs = -Z_constraint' * Mass_Matrix * B_constraint;
        force_reduced_rhs = Z_constraint' * Jv_end';
        reduced_rhs = [bias_reduced_rhs, input_reduced_rhs, force_reduced_rhs];
        % 直接构造 2x2 解析逆，防止符号求解器主动展开并长时间化简。
        det_M_reduced = M_reduced(1,1) * M_reduced(2,2) - ...
                        M_reduced(1,2) * M_reduced(2,1);
        inv_M_reduced = [ M_reduced(2,2), -M_reduced(1,2);
                         -M_reduced(2,1),  M_reduced(1,1)] / det_M_reduced;
        reduced_solution = inv_M_reduced * reduced_rhs;
    
        ddq_bias = -B_constraint * dotJa_dotQ + ...
                   Z_constraint * reduced_solution(:,1);
        ddq_input = B_constraint + Z_constraint * reduced_solution(:,2:3);
        ddq_force = Z_constraint * reduced_solution(:,4:6);
    elseif N==1
        % 单节时没有内部零空间自由度，约束直接确定广义加速度。
        ddq_bias = -B_constraint * dotJa_dotQ;
        ddq_input = B_constraint;
        ddq_force = sym(zeros(2,3));
    end
    
    X = [q; dotQ];
    Fx = [dotQ; ddq_bias];
    Gx = [zeros(2*N,2); ddq_input];
    Hx = [zeros(2*N,3); ddq_force];
    stateDot = Fx + Gx * u_acc + Hx * F_ext;
end
disp('模型推导完成！');
