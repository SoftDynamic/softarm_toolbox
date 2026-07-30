function [ddq,tension,isPullOnlyFeasible,reciprocalCondition] = softarm_actuator_acceleration(q,dq,coordinateAcceleration,tauExternal,w,p)
%SOFTARM_ACTUATOR_ACCELERATION Enforce independent actuator-coordinate accelerations.
%#codegen
assert(numel(coordinateAcceleration)==3);
coordinateAcceleration=coordinateAcceleration(:);
M=softarm_mass(q,p); h=softarm_bias(q,dq,p); J=softarm_end_jacobian(q,p);
Ja=softarm_actuator_jacobian(q,p);
jdotdq=softarm_actuator_velocity_bias(q,dq,p);
S=Ja*(M\Ja.'); reciprocalCondition=rcond(S);
assert(isfinite(reciprocalCondition) && reciprocalCondition>sqrt(eps),'softarm:SingularActuatorConstraint','Actuator acceleration constraints are singular.');
n=numel(q); m=size(Ja,1); solution=[M,Ja.';Ja,zeros(m)]\[tauExternal+J.'*w-h;coordinateAcceleration-jdotdq];
ddq=solution(1:n); tension=solution(n+1:n+m);
isPullOnlyFeasible=all(tension([1 2 3])>=-sqrt(eps));
end
