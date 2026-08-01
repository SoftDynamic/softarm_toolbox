function [poseFrames,baseFrames,centerlineFrames] = ...
    reconstructBackbonePoses(qFrames,plant,samplesPerSegment)
%RECONSTRUCTBACKBONEPOSES Reconstruct the arm root and section-end poses.
%   Each row of QFRAMES is one generalized-coordinate sample. POSEFRAMES
%   stacks the arm mounting root and section-end transforms. BASEFRAMES
%   stacks the vehicle-body transform. CENTERLINEFRAMES stacks 3-D points
%   beginning at the arm root and continuing section-by-section.
arguments
    qFrames (:,:) double
    plant (1,1) struct
    samplesPerSegment (1,1) double {mustBeInteger, ...
        mustBeGreaterThanOrEqual(samplesPerSegment,2), ...
        mustBeLessThanOrEqual(samplesPerSegment,128)} = 16
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
baseFrames = zeros(frameCount,16);
hasCenterline = isfield(plant,"centerline") && ~isempty(plant.centerline);
if hasCenterline
    sampleCoordinates = (1:samplesPerSegment)/samplesPerSegment;
    centerlinePointCount = 1+segmentCount*samplesPerSegment;
else
    warning("softarm:MissingCenterline", ...
        "Bundle has no model-derived centerline; using the key-node polyline.");
    sampleCoordinates = [];
    centerlinePointCount = segmentCount+1;
end
centerlineFrames = zeros(frameCount,3*centerlinePointCount);
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
    baseFrames(frame,:) = baseTransform(:).';
    sectionTransforms = plant.kinematics(q,plant.parameters);
    if size(sectionTransforms,1) ~= 4 || size(sectionTransforms,2) ~= 4 || ...
            size(sectionTransforms,3) ~= segmentCount
        error("softarm:InvalidKinematics", ...
            "Kinematics output must contain %d 4-by-4 section transforms.", ...
            segmentCount);
    end
    rootTransform = baseTransform*mountTransform;
    transforms = cat(3,rootTransform,sectionTransforms);
    poseFrames(frame,:) = transforms(:).';
    if hasCenterline
        sampledPositions = plant.centerline(q,plant.parameters,sampleCoordinates);
        if ~isequal(size(sampledPositions),[3,segmentCount*samplesPerSegment])
            error("softarm:InvalidCenterline", ...
                "Centerline output dimensions do not match the model and sample count.");
        end
        centerline = [rootTransform(1:3,4),sampledPositions];
    else
        centerline = reshape(transforms(1:3,4,:),3,segmentCount+1);
    end
    centerlineFrames(frame,:) = centerline(:).';
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
