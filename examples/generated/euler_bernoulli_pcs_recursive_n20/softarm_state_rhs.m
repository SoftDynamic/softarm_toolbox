function dx=softarm_state_rhs(x,tauArm,wVehicle,wTip,p)
%#codegen
n=numel(x)/2; q=x(1:n); dq=x(n+1:end); dx=[dq;softarm_forward_dynamics(q,dq,tauArm,wVehicle,wTip,p)];
end
