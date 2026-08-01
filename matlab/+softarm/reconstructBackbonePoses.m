function poseFrames = reconstructBackbonePoses(qFrames,plant)
%RECONSTRUCTBACKBONEPOSES Reconstruct the arm root and section-end poses.
%   POSEFRAMES = RECONSTRUCTBACKBONEPOSES(QFRAMES,PLANT) accepts one
%   generalized-coordinate sample per row. Each output row stacks the
%   world transform of the arm mounting root followed by every section-end
%   world transform, using MATLAB column-major order for each 4-by-4 pose.
arguments
    qFrames (:,:) double
    plant (1,1) struct
end

if ~all(isfield(plant,{'nq','nbase','parameters','manifest','kinematics'}))
    error("softarm:InvalidPlant","Plant metadata or function handles are incomplete.");
end
if size(qFrames,2) ~= plant.nq
    error("softarm:InvalidQFrames", ...
        "The q log width (%d) does not match bundle nq (%d).", ...
        size(qFrames,2),plant.nq);
end
if ~isfield(plant.manifest,"model") || ...
        ~all(isfield(plant.manifest.model,{'base_mode','segments', ...
        'mount_xyz','mount_rpy'}))
    error("softarm:MissingBaseMetadata", ...
        "Bundle manifest lacks base mount metadata; regenerate the bundle.");
end

model = plant.manifest.model;
mountTransform = transformRpy(model.mount_xyz,model.mount_rpy);
segmentCount = double(model.segments);
frameCount = size(qFrames,1);
poseFrames = zeros(frameCount,16*(segmentCount+1));
for frame = 1:frameCount
    q = qFrames(frame,:).';
    switch string(model.base_mode)
        case "fixed"
            baseTransform = eye(4);
        case "floating_rpy"
            if plant.nbase ~= 6
                error("softarm:InvalidPlant", ...
                    "A floating_rpy plant must have six base coordinates.");
            end
            baseTransform = transformRpy(q(1:3),q(4:6));
        otherwise
            error("softarm:UnsupportedBaseMode", ...
                "Unsupported base mode '%s'.",string(model.base_mode));
    end
    sectionTransforms = plant.kinematics(q,plant.parameters);
    if size(sectionTransforms,1) ~= 4 || size(sectionTransforms,2) ~= 4 || ...
            size(sectionTransforms,3) ~= segmentCount
        error("softarm:InvalidKinematics", ...
            "Kinematics output must contain %d 4-by-4 section transforms.", ...
            segmentCount);
    end
    transforms = cat(3,baseTransform*mountTransform,sectionTransforms);
    poseFrames(frame,:) = transforms(:).';
end
end

function transform = transformRpy(position,rpy)
position = reshape(double(position),3,1);
rpy = reshape(double(rpy),3,1);
roll = rpy(1); pitch = rpy(2); yaw = rpy(3);
cx = cos(roll); sx = sin(roll);
cy = cos(pitch); sy = sin(pitch);
cz = cos(yaw); sz = sin(yaw);
rotation = [ ...
    cz*cy, cz*sy*sx-sz*cx, cz*sy*cx+sz*sx; ...
    sz*cy, sz*sy*sx+cz*cx, sz*sy*cx-cz*sx; ...
    -sy, cy*sx, cy*cx];
transform = [rotation,position;0,0,0,1];
end
