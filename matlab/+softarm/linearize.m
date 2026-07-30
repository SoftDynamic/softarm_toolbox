function [A,Barm,Bvehicle,Btip] = linearize(plant, x, tauArm, wVehicle, wTip, p, step)
%LINEARIZE Central-difference linearization of the generated state RHS.
arguments
    plant (1,1) struct
    x (:,1) double
    tauArm (:,1) double
    wVehicle (6,1) double = zeros(6,1)
    wTip (6,1) double = zeros(6,1)
    p (:,1) double = plant.parameters
    step (1,1) double {mustBePositive} = eps^(1/3)
end
nx = numel(x);
nu = numel(tauArm);
A = zeros(nx, nx);
Barm = zeros(nx, nu);
Bvehicle = zeros(nx, 6);
Btip = zeros(nx, 6);
for column = 1:nx
    delta = step * max(1, abs(x(column)));
    dx = zeros(nx,1); dx(column) = delta;
    A(:,column) = (plant.stateRhs(x+dx,tauArm,wVehicle,wTip,p) - ...
                   plant.stateRhs(x-dx,tauArm,wVehicle,wTip,p))/(2*delta);
end
for column = 1:nu
    delta = step * max(1, abs(tauArm(column)));
    du = zeros(nu,1); du(column) = delta;
    Barm(:,column) = (plant.stateRhs(x,tauArm+du,wVehicle,wTip,p) - ...
                      plant.stateRhs(x,tauArm-du,wVehicle,wTip,p))/(2*delta);
end
for column = 1:6
    delta = step * max(1, abs(wVehicle(column)));
    dw = zeros(6,1); dw(column) = delta;
    Bvehicle(:,column) = (plant.stateRhs(x,tauArm,wVehicle+dw,wTip,p) - ...
                          plant.stateRhs(x,tauArm,wVehicle-dw,wTip,p))/(2*delta);
    delta = step * max(1, abs(wTip(column)));
    dw = zeros(6,1); dw(column) = delta;
    Btip(:,column) = (plant.stateRhs(x,tauArm,wVehicle,wTip+dw,p) - ...
                      plant.stateRhs(x,tauArm,wVehicle,wTip-dw,p))/(2*delta);
end
end
