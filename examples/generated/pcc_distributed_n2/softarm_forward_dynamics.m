function ddq = softarm_forward_dynamics(q,dq,tauArm,wVehicle,wTip,p)
%#codegen
M=softarm_mass(q,p); h=softarm_bias(q,dq,p); Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p); ddq=M\(Q-h);
end
