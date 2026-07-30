function ddq = softarm_forward_dynamics(q,dq,tau,w,p)
%#codegen
M=softarm_mass(q,p); h=softarm_bias(q,dq,p); J=softarm_end_jacobian(q,p); ddq=M\(tau+J.'*w-h);
end
