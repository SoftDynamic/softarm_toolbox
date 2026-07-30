function [ddq,tension,isPullOnlyFeasible,reciprocalCondition] = softarm_actuator_acceleration(q,dq,coordinateAcceleration,tauArmExternal,wVehicle,wTip,p)
%SOFTARM_ACTUATOR_ACCELERATION Enforce independent actuator-coordinate accelerations.
%#codegen
assert(numel(coordinateAcceleration)==3);
coordinateAcceleration=coordinateAcceleration(:);
M=softarm_mass(q,p); h=softarm_bias(q,dq,p);
Ja=softarm_actuator_jacobian(q,p);
JaFull=[zeros(3,0),Ja];
jdotdq=softarm_actuator_velocity_bias(q,dq,p);
S=JaFull*(M\JaFull.'); reciprocalCondition=rcond(S);
assert(isfinite(reciprocalCondition) && reciprocalCondition>sqrt(eps),'softarm:SingularActuatorConstraint','Actuator acceleration constraints are singular.');
Q=softarm_applied_force(q,tauArmExternal,wVehicle,wTip,p);
n=numel(q); m=size(JaFull,1); solution=[M,JaFull.';JaFull,zeros(m)]\[Q-h;coordinateAcceleration-jdotdq];
ddq=solution(1:n); tension=solution(n+1:n+m);
isPullOnlyFeasible=all(tension([1 2 3])>=-sqrt(eps));
end
