function y = softarm_actuator_coordinates(q,p)
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
act_pair_x_s1_radius = p(24);
act_pair_x_s2_radius = p(25);
act_pair_y_s1_radius = p(26);
act_pair_y_s2_radius = p(27);
t0 = 1./s1_length;
t1 = 1./s2_length;
y = reshape([-3*act_pair_x_s1_radius.*ax1.*t0/2 - 3*act_pair_x_s2_radius.*ax2.*t1/2;-act_pair_y_s1_radius.*t0.*(2.9398464800886433e-17*ax1 + 1.5*ay1) - act_pair_y_s2_radius.*t1.*(2.9398464800886433e-17*ax2 + 1.5*ay2)],2,1);
end
