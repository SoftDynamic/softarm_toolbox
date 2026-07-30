function plant = configureBundle(owner, bundle, requirement)
%CONFIGUREBUNDLE Resolve a Plant mask bundle and initialize its dimensions.
arguments
    owner (1,1) string
    bundle (1,1) string
    requirement (1,1) string = "plant"
end
model = bdroot(owner);
modelFile = string(get_param(model,"FileName"));
root = string(fileparts(modelFile));
candidate = bundle;
if ~java.io.File(char(candidate)).isAbsolute()
    candidate = fullfile(root,candidate);
end
candidate = string(java.io.File(char(candidate)).getCanonicalPath());
assert(isfolder(candidate),"softarm:MissingBundle", ...
    "Bundle directory does not exist: %s",candidate);
plant = softarm.initModel(candidate);
if requirement == "force" || requirement == "strict"
    assert(~isempty(plant.actuation),"softarm:MissingActuation", ...
        "Selected bundle has no actuator mapping.");
end
if requirement == "strict"
    assert(isfield(plant.actuation,"acceleration"), ...
        "softarm:MissingAcceleration", ...
        "Selected bundle has no strict acceleration mapping.");
end
set_param(owner,"UserData",struct( ...
    "bundle",candidate,"coordinates",{plant.manifest.coordinates}), ...
    "UserDataPersistent","on");
end
