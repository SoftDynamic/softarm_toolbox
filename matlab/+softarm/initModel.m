function plant = initModel(bundleDir)
%INITMODEL Initialize base-workspace data used by softarm_plant.slx.
arguments
    bundleDir (1,1) string
end
plant = softarm.loadModel(bundleDir);
assignin("base", "softarm_bundle", char(plant.bundle));
assignin("base", "softarm_nq", plant.nq);
assignin("base", "softarm_nbase", plant.nbase);
assignin("base", "softarm_narm", plant.narm);
assignin("base", "softarm_p", plant.parameters);
if ~isempty(plant.actuation)
    assignin("base", "softarm_nu", plant.actuation.count);
    assignin("base", "softarm_u0", zeros(plant.actuation.count,1));
end
q0 = zeros(plant.nq,1);
armCoordinates = cellstr(plant.manifest.coordinates.arm);
for armIndex = 1:plant.narm
    token = regexp(armCoordinates{armIndex}, '^l(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        parameterName = "s" + token{1} + "_rest_length";
        match = find(strcmp({plant.manifest.parameters.name}, parameterName), 1);
        if ~isempty(match), q0(plant.nbase+armIndex) = plant.parameters(match); end
    end
end
assignin("base", "softarm_x0", [q0;zeros(plant.nq,1)]);
if ~isempty(plant.constraint)
    assignin("base", "softarm_nc", plant.constraint.count);
    assignin("base", "softarm_ac0", zeros(plant.constraint.count,1));
end
end
