root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_kirchhoff_pcs_three_tendon_n2"));

q = [0;0;0.45;0;0;0.50];
dq = zeros(plant.nq,1);
coordinateAcceleration = zeros(plant.actuation.count,1);
[ddq,tension,isPullOnlyFeasible,reciprocalCondition] = ...
    plant.actuation.acceleration(q,dq,coordinateAcceleration, ...
        zeros(plant.narm,1),zeros(6,1),zeros(6,1),plant.parameters)
