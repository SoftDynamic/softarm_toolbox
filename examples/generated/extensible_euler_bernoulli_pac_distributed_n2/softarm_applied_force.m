function Q = softarm_applied_force(q,tauArm,wVehicle,wTip,p)
%SOFTARM_APPLIED_FORCE Assemble arm, vehicle-body, and world-tip loads.
%#codegen
assert(numel(tauArm)==8); assert(numel(wVehicle)==6); assert(numel(wTip)==6);
Q=[zeros(0,1);tauArm(:)]+softarm_vehicle_wrench_map(q,p)*wVehicle(:)+softarm_end_jacobian(q,p).'*wTip(:);
end
