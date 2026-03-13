clear; close all;
addpath("functions");

% 常数参量
N = 2; % PCC段数
g = 9.81; % NED坐标系

M_i = 0.2;                              % 软臂每节质量 (kg)
L_i_0 = 0.5;                            % 软臂每节原长 (m)
I_i_local = diag([0.002, 0.002, 1e-4]); % 软臂每节惯性矩
R = 0.0275;                             % 软臂半径 (m)
Me = 0.1;                               % 末端负载 (kg)
Ie_local = diag([1e-4, 1e-4, 1e-4]);    % 负载惯性矩
k_theta = 1.2;                          % 弯曲刚度 (N*m/rad)
k_l = 10;                               % 伸缩刚度 (N/m)
d_theta = 0.5;                          % 弯曲阻尼 (N*m*s/rad)
d_l = 5.0;                              % 轴向阻尼 (N*s/m)
d_phi = 0.1;                            % 扭转/方位角阻尼 (N*m*s/rad)

% 符号变量
syms theta phi l dotTheta dotPhi dotL [1 N] real % 广义坐标
syms T1 T2 T3 real                    % 绳上拉力
syms fx fy fz mx my mz real           % 末端外力与力矩
F_ext = [fx; fy; fz; mx; my; mz];     % 外部扳手 (Wrench)

%% 1. 运动学分析

Jw = sym(zeros(3,3*N,N));
Jv = sym(zeros(3,3*N,N));
Jv_com = sym(zeros(3,3*N,N));
q = reshape([theta;phi;l],[],1);
dotQ = reshape([dotTheta;dotPhi;dotL],[],1);
H_0 = sym(zeros(4,4,N));                            % H^0_i
H_0_com = sym(zeros(4,4,N));
w_0 = sym(zeros(3,1,N));                            % w^{0}_{0,i}

L1 = 0;
L2 = 0;
L3 = 0;

for i=1:N
    H_i = pcc_homogeneous(theta(i),phi(i),l(i));    % H^{i-1}_i
    H_com_i_local = pcc_homogeneous(theta(i)/2, phi(i), l(i)/2); % 中心（视为质心）齐次变换

    if i==1
        H_0(:,:,i) = H_i;
        H_0_com(:,:,i) = H_com_i_local;
    else
        H_0(:,:,i) = H_0(:,:,i-1) * H_i;
        H_0_com(:,:,i) = H_0(:,:,i-1) * H_com_i_local;
    end

    W = (diff(pcc_rotmat(theta(i),phi(i)), theta(i)) * dotTheta(i) + diff(pcc_rotmat(theta(i),phi(i)),phi(i)) * dotPhi(i))*(pcc_rotmat(theta(i),phi(i))');
    W = simplify(W);
    w_i = vex3(W);                                  % w^{i-1}_{i_1,i}
    if i==1
        w_0(:,:,i) = w_i;
    else
        w_0(:,:,i) = w_0(:,:,i-1) + H_0(1:3,1:3,i-1) * w_i; % w^{0}_{0,i-1} + w^{0}_{i_1,i}
    end
    
    Jw(:,:,i) = jacobian(w_0(:,:,i), dotQ);

    Jv(:,:,i) = jacobian(H_0(1:3,4,i), q);
    Jv_com(:,:,i) = jacobian(H_0_com(1:3,4,i), q);

    L1 = L1 + l(i) - theta(i)*R*cos(phi(i));
    L2 = L2 + l(i) - theta(i)*R*cos(phi(i)+2*sym(pi)/3);
    L3 = L3 + l(i) - theta(i)*R*cos(phi(i)-2*sym(pi)/3);
end

Ja = jacobian([L1;L2;L3], q);

%% 2. 能量项计算
% 计算末端雅可比 (Body Jacobian 或 Space Jacobian)
Jv_end = Jv(:,:,N);
Jw_end = Jw(:,:,N);
J_full_end = [Jv_end; Jw_end];

% 修正总动能 (包含机械臂各段 + 末端负载)
K_total = 0;
for i = 1:N
    % 1. 平动动能 (使用质心雅可比)
    K_trans = 0.5 * M_i * (dotQ' * (Jv_com(:,:,i)' * Jv_com(:,:,i)) * dotQ);
    
    % 2. 转动动能
    R_i = H_0_com(1:3, 1:3, i); % 当前段质心处的旋转矩阵
    I_i_base = R_i * I_i_local * R_i'; 
    K_rot = 0.5 * (dotQ' * Jw(:,:,i)' * I_i_base * Jw(:,:,i) * dotQ);
    
    K_total = K_total + K_trans + K_rot;
end
R_end = H_0(1:3, 1:3, N);
Ie_base = R_end * Ie_local * R_end'; % 变换到基坐标系
K_load = 0.5 * Me * (dotQ' * (Jv_end' * Jv_end) * dotQ) + ...
         0.5 * (dotQ' * Jw_end' * Ie_base * Jw_end * dotQ);
K_sum = K_total + K_load;

% 修正总势能 (重力 + TPU 弹性能)
P_gravity = 0;
for i = 1:N
    P_gravity = P_gravity - M_i * g * H_0_com(3,4,i);
end
P_load = - Me * g * H_0(3,4,N);
P_elastic = 0;
for i = 1:N
    % 弹性势能模型：弯曲 theta 和 伸缩 l 的回复力
    P_elastic = P_elastic + 0.5 * k_theta * theta(i)^2 + ...
                            0.5 * k_l * (l(i) - L_i_0)^2;
end
P_sum = P_gravity + P_load + P_elastic;

%% 3. 动力学矩阵提取
disp('正在计算质量矩阵...');
Mass_Matrix = simplify(jacobian(jacobian(K_sum, dotQ)', dotQ));

disp('正在计算重力与弹性矢量...');
G_K_Matrix = simplify(jacobian(P_sum, q)');

% 阻尼项
D_Matrix = diag(repmat([d_theta, d_phi, d_l], 1, N)) * dotQ;

% 科里奥利项 C(q, dotQ)
disp('正在计算科里奥利项...');

% % 初始化 C 矩阵 (3N x 3N)
% C_mat = sym(zeros(size(q,1), size(q,1)));
% 
% % 使用 Christoffel 符号计算: C_kj = sum_i( 0.5 * (dM_kj/dqi + dM_ki/dqj - dM_ij/dqk) * dotQi )
% for k = 1:size(q,1)
%     for j = 1:size(q,1)
%         for i = 1:size(q,1)
%             term = 0.5 * (diff(Mass_Matrix(k,j), q(i)) + ...
%                           diff(Mass_Matrix(k,i), q(j)) - ...
%                           diff(Mass_Matrix(i,j), q(k)));
%             C_mat(k,j) = C_mat(k,j) + term * dotQ(i);
%         end
%     end
% end
% 
% % 计算 C * dotQ 向量
% C_vector = simplify(C_mat * dotQ);

% 拉格朗日恒等式法
% C(q,dq)dq = dM*dq-\partial(1/2 dq' M dq)/\partial q
G_kinetic = jacobian(K_sum, q)'; % d/dq (1/2 * dq' * M * dq)

dotM_dq = sym(zeros(size(q,1), 1));
for i = 1:size(q,1)
    dotM_dq = dotM_dq + diff(Mass_Matrix, q(i)) * dotQ(i) * dotQ;
end

C_vector = simplify(dotM_dq - G_kinetic);

% 广义力映射
% Ja 是绳索拉力 T 到关节空间的映射 (3x3N)
% J_full_end 是末端外力到关节空间的映射 (6x3N)
disp('正在计算广义力...');
Tau_rope = Ja' * [T1; T2; T3];
Tau_ext = J_full_end' * F_ext;

%% 4. 导出为 MATLAB Function (用于 Simulink)
disp('正在导出函数文件...');

% 合力
% Total_Bias = C(q,dq)*dq + G_K(q) + D*dq - Tau_rope - Tau_ext
Total_Bias = C_vector + G_K_Matrix + D_Matrix - Tau_rope - Tau_ext;
Total_Bias_wo_Cori = G_K_Matrix + D_Matrix - Tau_rope - Tau_ext;

if not(isfolder('output'))
    mkdir 'output';
end

% 导出质量矩阵
matlabFunction(Mass_Matrix, 'File', 'output/get_MassMatrix', 'Vars', {q});

% 导出偏移力项 (包含科里奥利力、重力、弹性、阻尼及外力)
matlabFunction(Total_Bias, 'File', 'output/get_BiasForce', ...
    'Vars', {q, dotQ, T1, T2, T3, fx, fy, fz, mx, my, mz});

matlabFunction(Total_Bias_wo_Cori, 'File', 'output/get_BiasForceWoCori', ...
    'Vars', {q, dotQ, T1, T2, T3, fx, fy, fz, mx, my, mz});

disp('导出完成！');
