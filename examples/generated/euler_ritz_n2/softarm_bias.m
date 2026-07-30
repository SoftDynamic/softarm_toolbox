function h = softarm_bias(q,dq,p)
% Generated from the SymPy model. Do not edit.
%#codegen
ax1 = q(1);
ay1 = q(2);
ax2 = q(3);
ay2 = q(4);
dax1 = dq(1);
day1 = dq(2);
dax2 = dq(3);
day2 = dq(4);
s1_length = p(1);
s2_length = p(2);
s1_mass = p(3);
s2_mass = p(4);
s1_Ixx = p(5);
s2_Ixx = p(6);
s1_Iyy = p(7);
s2_Iyy = p(8);
s1_Izz = p(9);
s2_Izz = p(10);
s1_EI_x = p(11);
s2_EI_x = p(12);
s1_EI_y = p(13);
s2_EI_y = p(14);
s1_d_ax = p(15);
s2_d_ax = p(16);
s1_d_ay = p(17);
s2_d_ay = p(18);
gravity = p(19);
tip_mass = p(20);
tip_Ixx = p(21);
tip_Iyy = p(22);
tip_Izz = p(23);
t0 = 3./s1_length.^3;
t1 = 3./s2_length.^3;
h = reshape([ax1.*s1_EI_y.*t0 + dax1.*s1_d_ax;ay1.*s1_EI_x.*t0 + day1.*s1_d_ay;ax2.*s2_EI_y.*t1 + dax2.*s2_d_ax;ay2.*s2_EI_x.*t1 + day2.*s2_d_ay],4,1);
end
