root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "pcc_three_tendon_extensible_n2"));

q = [0;0;0.45;0;0;0.50];
dq = zeros(plant.nq,1);
tension = softarm.packActuatorInputs(plant.actuation, ...
    struct("t1",2.0,"t2",1.5,"t3",1.0));
[tau,isFeasible] = plant.actuation.force(q,tension,plant.parameters)
ddq = plant.forwardDynamics(q,dq,tau,zeros(6,1),plant.parameters)
