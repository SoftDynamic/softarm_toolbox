function Ja = softarm_actuator_jacobian(q,p)
% Generated from the SymPy model. Do not edit.
%#codegen
base_x = q(1);
base_y = q(2);
base_z = q(3);
base_roll = q(4);
base_pitch = q(5);
base_yaw = q(6);
bx1 = q(7);
by1 = q(8);
l1 = q(9);
s1_rest_length = p(1);
s1_mass = p(2);
s1_Ixx = p(3);
s1_Iyy = p(4);
s1_Izz = p(5);
s1_k_bx = p(6);
s1_k_by = p(7);
s1_k_l = p(8);
s1_d_bx = p(9);
s1_d_by = p(10);
s1_d_l = p(11);
gravity = p(12);
tip_mass = p(13);
tip_Ixx = p(14);
tip_Iyy = p(15);
tip_Izz = p(16);
vehicle_mass = p(17);
vehicle_Ixx = p(18);
vehicle_Iyy = p(19);
vehicle_Izz = p(20);
act_t1_s1_radius = p(21);
act_t2_s1_radius = p(22);
act_t3_s1_radius = p(23);
Ja = reshape([-act_t1_s1_radius;0.49999999999999983*act_t2_s1_radius;0.50000000000000042*act_t3_s1_radius;0;-0.86602540378443874*act_t2_s1_radius;0.8660254037844384*act_t3_s1_radius;1;1;1],3,3);
end
