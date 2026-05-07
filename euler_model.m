clear; close all;
addpath("functions"); 

%% 1. 参数定义
N = 2;              % 段数
g = 9.81;           % 重力加速度
L0 = 1.0 / N;           % 每节固定长度 (m)
R = 0.03;         % 驱动绳分布半径 (m)
Me = 0.1;           % 末端质点质量 (kg)
M = 0.3;            % 软臂质量

% % 刚度与阻尼
EI = 1.2;           % 弯曲刚度 (N*m^2) -> k_theta = EI/L0
d_x = 0.05;
d_y = 0.05;
% syms Me EI d_x d_y real
Mi = M / N ;
% sym_constants = {M, Me, EI, d_x, d_y};

% 符号变量定义
x = sym('x', [1 N], 'real');
y = sym('y', [1 N], 'real');
dotX = sym('dotX', [1 N], 'real'); % 关节速度 dx/dt
dotY = sym('dotY', [1 N], 'real'); % 关节速度 dy/dt

u_acc = sym('u_acc', [2 1], 'real');
% 末端外力 (仅考虑三维力，忽略质点力矩)
syms fx fy fz real
F_ext = [fx; fy; fz];

% 广义坐标合并
q = reshape([x; y], [], 1);
dotQ = reshape([dotX; dotY], [], 1);

%% 2. 运动学分析 (PCC)
H_0 = sym(zeros(4,4,N));

for i = 1:N
    Hi = taylor_homogeneous(x(i), y(i), L0);
    if i == 1
        H_0(:,:,i) = Hi;
    else
        H_0(:,:,i) = H_0(:,:,i-1) * Hi;
    end
end

% 末端位置 p_end
p_end = H_0(1:3, 4, N);

% 末端平动雅可比
Jv_end = jacobian(p_end, q);

%% 3. 绳索驱动映射 (4根绳简化为2个独立对)
% 在 x, y 坐标下: Lx = R*x, Ly = R*y
Lx_total = 0; Ly_total = 0;
for i = 1:N
    Lx_total = Lx_total + R * x(i);
    Ly_total = Ly_total + R * y(i);
end
Ja = jacobian([Lx_total; Ly_total], q);
dotJa_dotQ = jacobian(Ja * dotQ, q) * dotQ;

%% 4. 动力学矩阵提取
% 4.1 动能
K_sum = 0.5 * Me * (dotQ' * (Jv_end' * Jv_end) * dotQ);
if N > 1
    for i = 1:N-1
        Jv_i = jacobian(H_0(1:3,4,i), q);
        K_sum = K_sum + 0.5 * Mi * (dotQ' * (Jv_i' * Jv_i) * dotQ);
    end
end

% 4.2 势能 (末端重力势能 + 欧拉-伯努利弯曲能，软臂重力势能忽略)
P_gravity = - Me * [0, 0, -g] * p_end; % 假设 z 为竖直向上则用 +g，这里依NED惯例
k_bending = EI / L0;
P_elastic = 0;
for i = 1:N
    P_elastic = P_elastic + 0.5 * k_bending * (x(i)^2 + y(i)^2);
end
P_sum = P_gravity + P_elastic;

% 4.3 质量矩阵
disp('正在计算质量矩阵...');
Mass_Matrix = simplify(jacobian(jacobian(K_sum, dotQ)', dotQ));

% 4.4 重力与弹性力矢量
disp('正在计算保守力 (G+K)...');
G_K_Matrix = simplify(jacobian(P_sum, q)');

% 4.5 阻尼
D_Matrix = diag(repmat([d_x, d_y], 1, N)) * dotQ;

% 4.6 科里奥利项 (Christoffel 符号)
disp('正在计算科里奥利项...');
C_mat = sym(zeros(size(q,1)));
for k = 1:size(q,1)
    for j = 1:size(q,1)
        for i = 1:size(q,1)
            term = 0.5 * (diff(Mass_Matrix(k,j), q(i)) + ...
                          diff(Mass_Matrix(k,i), q(j)) - ...
                          diff(Mass_Matrix(i,j), q(k)));
            C_mat(k,j) = C_mat(k,j) + term * dotQ(i);
        end
    end
end
C_vector = simplify(C_mat * dotQ);

% 4.7 动力学方程 \ddot{q} = M^{-1} * (Tau_rope + Tau_ext - C - G_K - D)

disp('正在收集动力学方程...');

n_vec = C_vector + G_K_Matrix + D_Matrix; % 非线性项合并为 n(q, \dot{q})
Tau_ext = Jv_end' * F_ext;
A_mat = Ja * (Mass_Matrix \ Ja');
Tau_rope = Ja' * (A_mat \ (u_acc - dotJa_dotQ + Ja * (Mass_Matrix \ n_vec) - Ja * (Mass_Matrix \ Tau_ext)));

%% 5. 状态空间模型

X = [q; dotQ];
dotX = [dotQ; Mass_Matrix \ (Tau_rope + Tau_ext - n_vec)];
dotX = simplify(dotX);

Gx = jacobian(dotX, u_acc);
Hx = jacobian(dotX, F_ext);
Fx = simplify(dotX - Gx * u_acc - Hx * F_ext);


% %% 6. 导出函数
% disp('正在导出函数...');
% if ~isfolder('output'), mkdir 'output'; end
% 
% matlabFunction(Mass_Matrix, 'File', 'output/get_MassMatrix_Euler', 'Vars', [{q},sym_constants]);
% matlabFunction(C_vector,    'File', 'output/get_Coriolis_Euler', 'Vars', [{q, dotQ},sym_constants]);
% matlabFunction(G_K_Matrix,  'File', 'output/get_GKForce_Euler', 'Vars', [{q},sym_constants]);
% matlabFunction(D_Matrix, 'File', 'output/get_DampForce_Euler', 'Vars', [{dotQ},sym_constants]);
% matlabFunction(Tau_rope + Tau_ext, 'File', 'output/get_ExtForce_Euler', 'Vars', [{q, dotQ, u_acc, fx, fy, fz},sym_constants]);
% matlabFunction(H_0, 'File', 'output/get_homogeneous_Euler', 'Vars', [{q},sym_constants]);
% 
disp('模型导出完成！');
