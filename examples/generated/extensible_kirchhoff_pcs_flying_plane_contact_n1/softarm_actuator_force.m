function [tau,isFeasible] = softarm_actuator_force(q,tension,p)
%SOFTARM_ACTUATOR_FORCE Map ordered channel tensions to generalized force.
%#codegen
assert(numel(tension)==3);
tension=tension(:);
isFeasible=all(tension([1 2 3])>=0);
assert(isFeasible,'softarm:NegativeUnilateralTension','Unilateral tendon tension must be nonnegative.');
Ja=softarm_actuator_jacobian(q,p);
tau=-Ja.'*tension;
end
