function Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p)
%#codegen
assert(numel(tauArm)==4); assert(numel(wVehicle)==6); assert(numel(wTip)==6);
Q=[zeros(6,1);tauArm(:)]+softarm_vehicle_wrench_map(q,p)*wVehicle(:)+softarm_end_jacobian(q,p).'*wTip(:);
end
