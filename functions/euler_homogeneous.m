function H = euler_homogeneous(x,y,l)
% x = theta*cos(phi) y = theta * sin(phi)


arguments (Input)
    x   (1,1) {mustBeReal}
    y     (1,1) {mustBeReal}
    l       (1,1) {mustBeReal}
end
arguments (Output)
    H       (4,4) {mustBeReal}
end
H = [
[(x^2*cos((x^2 + y^2)^(1/2)) + y^2)/(x^2 + y^2), (x*y*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2), (x*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2), -(l*x*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2)]
[(x*y*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2), (y^2*cos((x^2 + y^2)^(1/2)) + x^2)/(x^2 + y^2), (y*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2), -(l*y*(cos((x^2 + y^2)^(1/2)) - 1))/(x^2 + y^2)]
[ -(x*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2),  -(y*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2),                       cos((x^2 + y^2)^(1/2)),    (l*sin((x^2 + y^2)^(1/2)))/(x^2 + y^2)^(1/2)]
[                                             0,                                              0,                                            0,                                               1]
];
end