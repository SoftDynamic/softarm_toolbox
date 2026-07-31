root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_kirchhoff_pcs_signed_pair_n2"));

q = [0;0;0.45;0;0;0.50];
tension = softarm.packActuatorInputs(plant.actuation, ...
    struct("bend_x",-2.0,"bend_y",1.5));
[tau,isFeasible] = plant.actuation.force(q,tension,plant.parameters)
