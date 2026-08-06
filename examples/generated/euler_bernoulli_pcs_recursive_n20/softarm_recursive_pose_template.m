function H = softarm_recursive_pose_template(qLocal,pLocal,xi)
% Generated from the SymPy model. Do not edit.
%#codegen
local_q1 = qLocal(1);
local_q2 = qLocal(2);
local_p1 = pLocal(1);
local_p2 = pLocal(2);
local_p3 = pLocal(3);
local_p4 = pLocal(4);
local_p5 = pLocal(5);
local_p6 = pLocal(6);
local_p7 = pLocal(7);
local_p8 = pLocal(8);
local_p9 = pLocal(9);
t0 = local_p1.^2;
t1 = xi.^2;
t2 = 1./t0;
t3 = local_q1.^2;
t4 = local_q2.^2;
t5 = t0.*t1.*(t2.*t3 + t2.*t4);
t6 = softarm_cosc_sqrt(t5);
t7 = t1.*t3;
t8 = local_q1.*t1.*t6;
t9 = -local_q2.*t8;
t10 = xi.*softarm_sinc_sqrt(t5);
t11 = local_q1.*t10;
t12 = t1.*t4;
t13 = local_q2.*t10;
t14 = -t12 - t7;
H = reshape([-t6.*t7 + 1;t9;-t11;0;t9;-t12.*t6 + 1;-t13;0;t11;t13;t14.*t6 + 1;0;local_p1.*t8;local_p1.*local_q2.*t1.*t6;local_p1.*xi.*(t14.*softarm_sinc3_sqrt(t5) + 1);1],4,4);
end
