function R = pcc_rotmat(theta,phi)

arguments (Input)
    theta   (1,1) {mustBeReal}
    phi     (1,1) {mustBeReal}
end
arguments (Output)
    R       (3,3) {mustBeReal}
end

Rz_posPhi  = rot2dcm([0 0 1], phi);   % 绕  +Z 轴转  phi
Ry_theta   = rot2dcm([0 1 0], theta);   % 绕  +Y 轴转 theta
Rz_negPhi  = rot2dcm([0 0 1], -phi);    % 绕  –Z 轴转  phi

% 合成总旋转矩阵
R = Rz_posPhi * Ry_theta * Rz_negPhi;
end