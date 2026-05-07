function H = taylor_homogeneous(x,y,l)
% x = theta*cos(phi) y = theta * sin(phi)


arguments (Input)
    x   (1,1) {mustBeReal}
    y     (1,1) {mustBeReal}
    l       (1,1) {mustBeReal}
end
arguments (Output)
    H       (4,4) {mustBeReal}
end

theta_sq = x^2 + y^2;

f1 = 1 - theta_sq/6;% + theta_sq^2/120;           % Taylor of sin(theta)/theta
f2 = 1/2 - theta_sq/24;% + theta_sq^2/720;        % Taylor of (1-cos(theta))/theta^2


R = [ 1 - y^2*f2,   x*y*f2,       x*f1;
      x*y*f2,       1 - x^2*f2,   y*f1;
     -x*f1,        -y*f1,         1 - theta_sq*f2 ];
p = [ l*x*f2; 
      l*y*f2; 
      l*f1 ];

H = [R, p; 0 0 0 1];

% H = [
% [(x^2*cos((x^2 + y^2)^(1/2)) + y^2)/(x^2 + y^2), (x*y*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2), (x*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2), -(l*x*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2)]
% [(x*y*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2), (y^2*cos((x^2 + y^2)^(1/2)) + x^2)/(x^2 + y^2), (y*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2), -(l*y*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2)]
% [ -(x*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2),  -(y*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2),                       cos((x^2 + y^2)^(1/2)),    (l*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2)]
% [                                             0,                                              0,                                            0,                                               1]
% ];
end