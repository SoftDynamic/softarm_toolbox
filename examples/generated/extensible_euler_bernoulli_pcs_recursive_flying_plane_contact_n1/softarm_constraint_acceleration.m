function [ddq,reaction,isFeasible,reciprocalCondition] = softarm_constraint_acceleration(q,dq,tauArm,wVehicle,wTip,constraintAcceleration,p)
%#codegen
assert(numel(constraintAcceleration)==1);
assert(norm(p([30 31 32]))>sqrt(eps),'softarm:InvalidPlaneNormal','Plane normal must be nonzero.');
constraintAcceleration=constraintAcceleration(:);
M=softarm_mass(q,p); h=softarm_bias(q,dq,p);
Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p);
[phi,A,gamma,G]=softarm_constraint_terms(q,dq,p);
gains=softarm_constraint_stabilization(p); omega=gains(:,1); zeta=gains(:,2);
targetAcceleration=constraintAcceleration-gamma-2*zeta.*omega.*(A*dq)-(omega.^2).*phi;
n=numel(q); m=size(A,1); K=[M,-G;A,zeros(m)]; reciprocalCondition=rcond(K);
assert(isfinite(reciprocalCondition)&&reciprocalCondition>eps,'softarm:SingularConstraint','Constraint KKT system is singular or ill-conditioned.');
solution=K\[Q-h;targetAcceleration]; ddq=solution(1:n); reaction=solution(n+1:n+m);
isFeasible=all(reaction([1])>=-sqrt(eps));
end
