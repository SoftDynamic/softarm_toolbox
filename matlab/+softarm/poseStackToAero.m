function [translation,rpy] = poseStackToAero(backbonePoses)
%POSESTACKTOAERO Convert stacked world transforms to AERO actor transforms.
%   BACKBONEPOSES stores one 4-by-4 homogeneous transform per node using
%   MATLAB column-major order. TRANSLATION and RPY contain one row per node;
%   RPY uses roll-pitch-yaw angles whose rotation matrix is Rz*Ry*Rx.
%#codegen
arguments
    backbonePoses (:,1) double
end

elementCount = numel(backbonePoses);
assert(elementCount >= 16 && rem(elementCount,16) == 0, ...
    "softarm:InvalidBackbonePoses", ...
    "backbone_poses must contain a positive whole number of 4-by-4 transforms.");
nodeCount = elementCount/16;
transforms = reshape(backbonePoses,4,4,nodeCount);
translation = zeros(nodeCount,3);
rpy = zeros(nodeCount,3);

for node = 1:nodeCount
    transform = transforms(:,:,node);
    translation(node,:) = transform(1:3,4).';
    rotation = transform(1:3,1:3);
    horizontal = hypot(rotation(1,1),rotation(2,1));
    pitch = atan2(-rotation(3,1),horizontal);
    if horizontal > 1e-9
        roll = atan2(rotation(3,2),rotation(3,3));
        yaw = atan2(rotation(2,1),rotation(1,1));
    else
        % At gimbal lock, fix yaw at zero and retain an equivalent roll.
        roll = atan2(-rotation(2,3),rotation(2,2));
        yaw = 0;
    end
    rpy(node,:) = [roll,pitch,yaw];
end
end
