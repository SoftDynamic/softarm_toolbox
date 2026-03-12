function R = rot2dcm(axis, ang)
% 轴角 -> 旋转矩阵（Rodrigues）
    axis = axis/norm(axis);          % 单位化
    [x,y,z] = deal(axis(1),axis(2),axis(3));
    c = cos(ang); s = sin(ang); C = 1-c;
    R = [ x*x*C+c    x*y*C-z*s  x*z*C+y*s;
          y*x*C+z*s  y*y*C+c    y*z*C-x*s;
          z*x*C-y*s  z*y*C+x*s  z*z*C+c   ];
end