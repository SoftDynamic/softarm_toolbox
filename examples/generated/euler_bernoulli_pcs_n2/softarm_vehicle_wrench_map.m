function Bv = softarm_vehicle_wrench_map(q,p)
% Generated from the SymPy model. Do not edit.
%#codegen
bx1 = q(1);
by1 = q(2);
bx2 = q(3);
by2 = q(4);
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
s1_d_bx = p(15);
s2_d_bx = p(16);
s1_d_by = p(17);
s2_d_by = p(18);
gravity = p(19);
tip_mass = p(20);
tip_Ixx = p(21);
tip_Iyy = p(22);
tip_Izz = p(23);
Bv = reshape([0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0],4,6);
end
