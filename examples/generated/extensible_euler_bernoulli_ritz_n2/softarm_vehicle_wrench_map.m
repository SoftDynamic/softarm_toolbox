function Bv = softarm_vehicle_wrench_map(q,p)
% Generated from the SymPy model. Do not edit.
%#codegen
ax1 = q(1);
ay1 = q(2);
az1 = q(3);
ax2 = q(4);
ay2 = q(5);
az2 = q(6);
s1_rest_length = p(1);
s2_rest_length = p(2);
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
s1_EA = p(15);
s2_EA = p(16);
s1_d_ax = p(17);
s2_d_ax = p(18);
s1_d_ay = p(19);
s2_d_ay = p(20);
s1_d_az = p(21);
s2_d_az = p(22);
gravity = p(23);
tip_mass = p(24);
tip_Ixx = p(25);
tip_Iyy = p(26);
tip_Izz = p(27);
Bv = reshape([0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0],6,6);
end
