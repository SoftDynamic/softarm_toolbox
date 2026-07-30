function [A,B,E] = linearize(plant, x, tau, wrench, p, step)
%LINEARIZE Central-difference linearization of the generated state RHS.
arguments
    plant (1,1) struct
    x (:,1) double
    tau (:,1) double
    wrench (6,1) double = zeros(6,1)
    p (:,1) double = plant.parameters
    step (1,1) double {mustBePositive} = eps^(1/3)
end
nx = numel(x);
nu = numel(tau);
A = zeros(nx, nx);
B = zeros(nx, nu);
E = zeros(nx, 6);
for column = 1:nx
    delta = step * max(1, abs(x(column)));
    dx = zeros(nx,1); dx(column) = delta;
    A(:,column) = (plant.stateRhs(x+dx,tau,wrench,p) - ...
                   plant.stateRhs(x-dx,tau,wrench,p))/(2*delta);
end
for column = 1:nu
    delta = step * max(1, abs(tau(column)));
    du = zeros(nu,1); du(column) = delta;
    B(:,column) = (plant.stateRhs(x,tau+du,wrench,p) - ...
                   plant.stateRhs(x,tau-du,wrench,p))/(2*delta);
end
for column = 1:6
    delta = step * max(1, abs(wrench(column)));
    dw = zeros(6,1); dw(column) = delta;
    E(:,column) = (plant.stateRhs(x,tau,wrench+dw,p) - ...
                   plant.stateRhs(x,tau,wrench-dw,p))/(2*delta);
end
end
