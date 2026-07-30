function J = softarm_end_jacobian(q,p)
% Generated from the SymPy model. Do not edit.
%#codegen
ax1 = q(1);
ay1 = q(2);
ax2 = q(3);
ay2 = q(4);
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
t0 = 3./(2*s1_length);
t1 = s2_length.*t0 + 1;
t2 = 3./(2*s2_length);
J = reshape([t1;0;0;0;t0;0;0;t1;0;-t0;0;0;1;0;0;0;t2;0;0;1;0;-t2;0;0],6,4);
end
