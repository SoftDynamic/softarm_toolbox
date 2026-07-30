function plant = initModel(bundleDir)
%INITMODEL Initialize base-workspace data used by softarm_plant.slx.
arguments
    bundleDir (1,1) string
end
plant = softarm.loadModel(bundleDir);
assignin("base", "softarm_bundle", char(plant.bundle));
assignin("base", "softarm_nq", plant.nq);
assignin("base", "softarm_p", plant.parameters);
if ~isempty(plant.actuation)
    assignin("base", "softarm_nu", plant.actuation.count);
    assignin("base", "softarm_u0", zeros(plant.actuation.count,1));
end
q0 = zeros(plant.nq,1);
for coordinateIndex = 1:plant.nq
    token = regexp(plant.manifest.coordinates{coordinateIndex}, '^l(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        parameterName = "s" + token{1} + "_rest_length";
        match = find(strcmp({plant.manifest.parameters.name}, parameterName), 1);
        if ~isempty(match), q0(coordinateIndex) = plant.parameters(match); end
    end
end
assignin("base", "softarm_x0", [q0;zeros(plant.nq,1)]);
end
