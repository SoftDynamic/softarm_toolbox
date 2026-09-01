function y = softarm_legacy_affine_template(pLocal,xi)
% Generated from the SymPy model. Do not edit.
%#codegen
local_p1 = pLocal(1);
local_p2 = pLocal(2);
local_p3 = pLocal(3);
local_p4 = pLocal(4);
local_p5 = pLocal(5);
local_p6 = pLocal(6);
local_p7 = pLocal(7);
local_p8 = pLocal(8);
local_p9 = pLocal(9);
t0 = 3*xi.^2/2;
t1 = (-t0 + 3*xi)./local_p1;
t2 = -t1;
t3 = t0 - xi.^3/2;
y = reshape([1;0;0;0;0;1;0;0;0;0;1;0;0;0;local_p1.*xi;1;0;0;t2;0;0;0;0;0;t1;0;0;0;t3;0;0;0;0;0;0;0;0;0;t2;0;0;t1;0;0;0;t3;0;0],48,1);
end
