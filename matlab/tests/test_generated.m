function tests = test_generated
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
addpath(fullfile(root,"matlab","simulink"));
addpath(fullfile(root,"examples","simulink"));
testCase.TestData.root = root;
end

function testReferenceBundles(testCase)
names = [
    "extensible_euler_bernoulli_pcs_lumped_n2", ...
    "extensible_euler_bernoulli_pcs_distributed_n2", ...
    "euler_bernoulli_pac_distributed_n2", ...
    "extensible_euler_bernoulli_pac_distributed_n2", ...
    "euler_bernoulli_ritz_n2","euler_bernoulli_pcs_n2", ...
    "extensible_euler_bernoulli_ritz_n2", ...
    "cosserat_pcs_lumped_n1","cosserat_pcs_distributed_tendon_n1"
];
for name = names
    bundle = fullfile(testCase.TestData.root,"examples","generated",name);
    plant = softarm.loadModel(bundle);
    q = zeros(plant.nq,1);
    armCoordinates = cellstr(plant.manifest.coordinates.arm);
    for k = 1:plant.narm
        token = regexp(armCoordinates{k},'^l(\d+)$','tokens','once');
        if ~isempty(token)
            field = "s"+token{1}+"_rest_length";
            index = find(strcmp({plant.manifest.parameters.name},field),1);
            q(plant.nbase+k) = plant.parameters(index);
        end
    end
    M = plant.mass(q,plant.parameters);
    verifyTrue(testCase,all(isfinite(M(:))));
    verifyLessThan(testCase,norm(M-M','fro'),1e-8);
    verifyGreaterThan(testCase,min(eig((M+M')/2)),0);
    dx = plant.stateRhs([q;zeros(plant.nq,1)],zeros(plant.narm,1), ...
        zeros(6,1),zeros(6,1),plant.parameters);
    verifyTrue(testCase,all(isfinite(dx)));
end
end

function testCosseratReferenceAndTendon(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "cosserat_pcs_distributed_tendon_n1"));
q = zeros(plant.nq,1);
dq = zeros(plant.nq,1);
verifyEqual(testCase,cellstr(plant.manifest.coordinates.arm), ...
    {'kx1';'ky1';'kz1';'vx1';'vy1';'vz1'});
H = plant.kinematics(q,plant.parameters);
verifyEqual(testCase,H(:,:,1),[eye(3),[0;0;0.45];0,0,0,1],"AbsTol",1e-12);
lengths = plant.actuation.coordinates(q,plant.parameters);
verifyEqual(testCase,lengths,0.45*ones(3,1),"AbsTol",1e-12);
Ja = plant.actuation.jacobian(q,plant.parameters);
verifyEqual(testCase,rank(Ja),3);
tension = [1.0;1.5;2.0];
[tau,isFeasible] = plant.actuation.force(q,tension,plant.parameters);
verifyTrue(testCase,isFeasible);
verifyEqual(testCase,tau,-Ja.'*tension,"AbsTol",1e-12);
ddq = plant.forwardDynamics(q,dq,tau,zeros(6,1),zeros(6,1),plant.parameters);
verifyTrue(testCase,all(isfinite(ddq)));
end

function testLinearization(testCase)
bundle = fullfile(testCase.TestData.root,"examples","generated","euler_bernoulli_ritz_n2");
plant = softarm.loadModel(bundle);
x = zeros(2*plant.nq,1); tau = zeros(plant.narm,1);
[A,Barm,Bvehicle,Btip] = softarm.linearize(plant,x,tau);
verifySize(testCase,A,[2*plant.nq,2*plant.nq]);
verifySize(testCase,Barm,[2*plant.nq,plant.narm]);
verifySize(testCase,Bvehicle,[2*plant.nq,6]);
verifySize(testCase,Btip,[2*plant.nq,6]);
verifyTrue(testCase,all(isfinite([A(:);Barm(:);Bvehicle(:);Btip(:)])));
end

function testPoseStackToAero(testCase)
angles = [0 0 0; 0.25 -0.35 0.45; -0.4 pi/2 0.2];
translations = [0 0 0; 1 -2 3; -0.2 0.3 0.4];
transforms = repmat(eye(4),1,1,size(angles,1));
for node = 1:size(angles,1)
    transforms(1:3,1:3,node) = rotationZYX(angles(node,:));
    transforms(1:3,4,node) = translations(node,:).';
end
[actualTranslations,actualAngles] = softarm.poseStackToAero(transforms(:));
verifyEqual(testCase,actualTranslations,translations,"AbsTol",1e-12);
verifyTrue(testCase,all(isfinite(actualAngles(:))));
for node = 1:size(angles,1)
    verifyEqual(testCase,rotationZYX(actualAngles(node,:)), ...
        transforms(1:3,1:3,node),"AbsTol",1e-10);
end
verifyError(testCase,@() softarm.poseStackToAero(zeros(15,1)), ...
    "softarm:InvalidBackbonePoses");
end

function testReconstructBackbonePoses(testCase)
mountXYZ = [0.12 -0.23 0.34];
mountRPY = [0.41 -0.32 0.13];
floating = syntheticPlaybackPlant("floating_rpy",mountXYZ,mountRPY,16,6);
qFrames = zeros(2,16);
qFrames(1,1:6) = [1.1 -2.2 3.3 0.21 -0.31 0.42];
qFrames(2,1:6) = [-0.4 0.5 -0.6 -0.12 0.23 -0.34];
[poses,baseFrames,centerlineFrames] = ...
    softarm.reconstructBackbonePoses(qFrames,floating,4);
verifySize(testCase,poses,[2 32]);
verifySize(testCase,baseFrames,[2 16]);
verifySize(testCase,centerlineFrames,[2 15]);
for frame = 1:2
    transforms = reshape(poses(frame,:),4,4,2);
    actualBase = reshape(baseFrames(frame,:),4,4);
    expectedBase = transformZYX(qFrames(frame,1:3),qFrames(frame,4:6));
    expectedRoot = transformZYX(qFrames(frame,1:3),qFrames(frame,4:6))* ...
        transformZYX(mountXYZ,mountRPY);
    verifyEqual(testCase,actualBase,expectedBase,"AbsTol",1e-12);
    verifyEqual(testCase,transforms(:,:,1),expectedRoot,"AbsTol",1e-12);
    verifyEqual(testCase,transforms(:,:,2),eye(4),"AbsTol",1e-12);
    centerline = reshape(centerlineFrames(frame,:),3,5);
    verifyEqual(testCase,centerline(:,1),expectedRoot(1:3,4),"AbsTol",1e-12);
end

fixed = syntheticPlaybackPlant("fixed",mountXYZ,mountRPY,2,0);
fixedPoses = softarm.reconstructBackbonePoses([0.1 0.2],fixed);
fixedTransforms = reshape(fixedPoses,4,4,2);
verifyEqual(testCase,fixedTransforms(:,:,1), ...
    transformZYX(mountXYZ,mountRPY),"AbsTol",1e-12);

fallback = floating;
fallback.centerline = [];
lastwarn("");
[~,~,fallbackCenterline] = softarm.reconstructBackbonePoses(qFrames,fallback,4);
[~,warningId] = lastwarn;
verifyEqual(testCase,string(warningId),"softarm:MissingCenterline");
verifySize(testCase,fallbackCenterline,[2 6]);

verifyError(testCase,@() softarm.reconstructBackbonePoses(zeros(1,15),floating), ...
    "softarm:InvalidQFrames");
legacy = floating;
legacy.manifest.model = rmfield(legacy.manifest.model,"mount_xyz");
verifyError(testCase,@() softarm.reconstructBackbonePoses(qFrames,legacy), ...
    "softarm:MissingBaseMetadata");
end

function testFloatingBundleBackboneRoot(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_flying_plane_contact_n1"));
q = zeros(1,plant.nq);
q(1:6) = [1.1 -2.2 3.3 0.21 -0.31 0.42];
lengthIndex = find(strcmp({plant.manifest.parameters.name},"s1_rest_length"),1);
q(plant.nbase+3) = plant.parameters(lengthIndex);
[poses,baseFrames,centerlineFrames] = ...
    softarm.reconstructBackbonePoses(q,plant);
transforms = reshape(poses,4,4,plant.manifest.model.segments+1);
baseTransform = reshape(baseFrames,4,4);
verifyEqual(testCase,transforms(:,:,1), ...
    transformZYX(q(1:3),q(4:6)),"AbsTol",1e-12);
verifyEqual(testCase,baseTransform, ...
    transformZYX(q(1:3),q(4:6)),"AbsTol",1e-12);
verifyEqual(testCase,transforms(:,:,2:end), ...
    plant.kinematics(q.',plant.parameters),"AbsTol",1e-12);
verifySize(testCase,centerlineFrames,[1 51]);
centerline = reshape(centerlineFrames,3,17);
verifyEqual(testCase,centerline(:,end),transforms(1:3,4,end),"AbsTol",1e-12);
end

function testGeneratedModelCenterlines(testCase)
root = testCase.TestData.root;
names = ["extensible_euler_bernoulli_pcs_lumped_n2", ...
    "euler_bernoulli_ritz_n2","euler_bernoulli_pac_distributed_n2", ...
    "extensible_euler_bernoulli_pac_distributed_n2","cosserat_pcs_lumped_n1"];
sampleCoordinates = [0.25 0.5 0.75 1.0];
for name = names
    plant = softarm.loadModel(fullfile(root,"examples","generated",name));
    verifyTrue(testCase,logical(plant.manifest.model.centerline));
    verifyNotEmpty(testCase,plant.centerline);
    q = nominalConfiguration(plant);
    q(1) = q(1)+0.12;
    points = plant.centerline(q,plant.parameters,sampleCoordinates);
    verifySize(testCase,points,[3 plant.manifest.model.segments*4]);
    verifyTrue(testCase,all(isfinite(points(:))));
    endpoints = plant.kinematics(q,plant.parameters);
    verifyEqual(testCase,points(:,4:4:end), ...
        reshape(endpoints(1:3,4,:),3,[]),"AbsTol",1e-10);
end

plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_lumped_n2"));
q = nominalConfiguration(plant);
q(1) = 0.35;
points = plant.centerline(q,plant.parameters,(1:16)/16);
firstSection = points(:,1:16);
chords = diff([zeros(3,1),firstSection],1,2);
verifyGreaterThan(testCase,norm(cross(chords(:,1),chords(:,end))),1e-5);
verifyError(testCase,@() plant.centerline(q,plant.parameters,[-0.1 1]), ...
    "softarm:InvalidMaterialCoordinate");
end

function testPACReferenceGeometry(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "euler_bernoulli_pac_distributed_n2"));
q = nominalConfiguration(plant);
q(1:3) = [0.4;-0.2;0.3];
xi = [0.25 0.5 1.0];
points = plant.centerline(q,plant.parameters,xi);
lengthIndex = find(strcmp({plant.manifest.parameters.name},"s1_length"),1);
lengthValue = plant.parameters(lengthIndex);
sineIntegral = arrayfun(@(x) integral( ...
    @(v) sin(q(1)*v+0.5*q(2)*v.^2),0,x),xi);
cosineIntegral = arrayfun(@(x) integral( ...
    @(v) cos(q(1)*v+0.5*q(2)*v.^2),0,x),xi);
expected = lengthValue*[cos(q(3))*sineIntegral; ...
    sin(q(3))*sineIntegral;cosineIntegral];
verifyEqual(testCase,points(:,1:3),expected,"AbsTol",2e-12);

transforms = plant.kinematics(q,plant.parameters);
alpha = q(1)+q(2)/2;
rz = [cos(q(3)) -sin(q(3)) 0;sin(q(3)) cos(q(3)) 0;0 0 1];
ry = [cos(alpha) 0 sin(alpha);0 1 0;-sin(alpha) 0 cos(alpha)];
verifyEqual(testCase,transforms(1:3,1:3,1),rz*ry,"AbsTol",2e-12);
end

function testPlaybackRejectsPoseStack(testCase)
bundle = fullfile(testCase.TestData.root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_flying_plane_contact_n1");
legacyPoseLog = [0 zeros(1,16)];
verifyError(testCase,@() softarm.playbackBackbonePoses( ...
    legacyPoseLog,Bundle=string(bundle)),"softarm:InvalidQFrames");
end

function testPlaybackGraphicsAndSampling(testCase)
bundle = fullfile(testCase.TestData.root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_flying_plane_contact_n1");
plant = softarm.loadModel(bundle);
q = zeros(1,plant.nq);
q(1:6) = [10 20 30 0.2 -0.3 0.4];
q(7) = 0.35;
q(9) = plant.parameters(1);
viewer = softarm.playbackBackbonePoses([0 q],Bundle=string(bundle), ...
    AxisLength=0.05,SamplesPerSegment=8);
cleanup = onCleanup(@() closeIfValid(viewer));

backbone = findobj(viewer,"Tag","softarm_backbone");
nodes = findobj(viewer,"Tag","softarm_backbone_nodes");
nodeXAxis = findobj(viewer,"Tag","softarm_x_axes");
baseXAxis = findobj(viewer,"Tag","softarm_base_x_axis");
verifyEqual(testCase,numel(backbone.XData),9);
verifyEqual(testCase,numel(nodes.XData),2);
verifyEqual(testCase,norm([baseXAxis.UData baseXAxis.VData baseXAxis.WData]), ...
    0.1,"AbsTol",1e-12);
verifyEqual(testCase,norm([nodeXAxis.UData(1) nodeXAxis.VData(1) ...
    nodeXAxis.WData(1)]),0.05,"AbsTol",1e-12);
verifyGreaterThanOrEqual(testCase,min(backbone.XData),backbone.Parent.XLim(1));
verifyLessThanOrEqual(testCase,max(backbone.XData),backbone.Parent.XLim(2));
verifyGreaterThanOrEqual(testCase,min(backbone.YData),backbone.Parent.YLim(1));
verifyLessThanOrEqual(testCase,max(backbone.YData),backbone.Parent.YLim(2));
verifyGreaterThanOrEqual(testCase,min(backbone.ZData),backbone.Parent.ZLim(1));
verifyLessThanOrEqual(testCase,max(backbone.ZData),backbone.Parent.ZLim(2));
verifyError(testCase,@() softarm.playbackBackbonePoses( ...
    [0 q],Bundle=string(bundle),SamplesPerSegment=1), ...
    "MATLAB:validators:mustBeGreaterThanOrEqual");
clear cleanup
end

function testGeneratedPoseStacks(testCase)
root = testCase.TestData.root;
names = ["extensible_euler_bernoulli_pcs_three_tendon_n2","euler_bernoulli_ritz_n2"];
for name = names
    plant = softarm.loadModel(fullfile(root,"examples","generated",name));
    q = zeros(plant.nq,1);
    armCoordinates = cellstr(plant.manifest.coordinates.arm);
    for armIndex = 1:plant.narm
        token = regexp(armCoordinates{armIndex},'^l(\d+)$','tokens','once');
        if ~isempty(token)
            parameterName = "s"+token{1}+"_rest_length";
            parameterIndex = find(strcmp({plant.manifest.parameters.name},parameterName),1);
            q(plant.nbase+armIndex) = plant.parameters(parameterIndex);
        end
    end
    transforms = plant.kinematics(q,plant.parameters);
    [translation,rpy] = softarm.poseStackToAero(transforms(:));
    verifySize(testCase,translation,[plant.manifest.model.segments 3]);
    verifySize(testCase,rpy,[plant.manifest.model.segments 3]);
    verifyEqual(testCase,translation,permute(transforms(1:3,4,:),[3 1 2]), ...
        "AbsTol",1e-12);
    verifyTrue(testCase,all(isfinite(rpy(:))));
end
end

function testThreeTendonForceAndAxialMode(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_three_tendon_n2"));
verifyEqual(testCase,plant.actuation.count,3);
verifyEqual(testCase,cellstr(plant.actuation.names),{'t1','t2','t3'});
q = [0.01;-0.02;0.46;0.03;0.01;0.51];
dq = zeros(plant.nq,1);
tension = softarm.packActuatorInputs(plant.actuation, ...
    struct("t1",2.0,"t2",1.5,"t3",1.0));
w = zeros(6,1);
[tau,isFeasible] = plant.actuation.force(q,tension,plant.parameters);
Ja = plant.actuation.jacobian(q,plant.parameters);
verifyTrue(testCase,isFeasible);
verifyEqual(testCase,tau,-Ja.'*tension,"AbsTol",1e-12);
verifyError(testCase,@() plant.actuation.force(q,[-1;0;0],plant.parameters), ...
    "softarm:NegativeUnilateralTension");
verifyEqual(testCase,Ja(:,[3,6]),ones(3,2),"AbsTol",1e-12);
verifyEqual(testCase,sum(Ja(:,[1,2,4,5]),1),zeros(1,4),"AbsTol",1e-12);
ddq = plant.forwardDynamics(q,dq,tau,zeros(6,1),w,plant.parameters);
verifyTrue(testCase,all(isfinite(ddq)));
deltaQ = [0.003;-0.002;0.001;0.004;-0.001;0.002];
verifyEqual(testCase,tau.'*deltaQ,-tension.'*(Ja*deltaQ),"AbsTol",1e-12);
end

function testThreeTendonAccelerationReference(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_three_tendon_n2"));
q = [0.01;-0.02;0.46;0.03;0.01;0.51];
dq = [0.02;0.01;0.001;-0.01;0.03;-0.002];
command = [0.01;-0.01;0.005];
w = zeros(6,1);
[ddq,tension,isPullOnlyFeasible,reciprocalCondition] = ...
    plant.actuation.acceleration(q,dq,command,zeros(plant.narm,1), ...
        zeros(6,1),w,plant.parameters);
Ja = plant.actuation.jacobian(q,plant.parameters);
jdotdq = plant.actuation.velocityBias(q,dq,plant.parameters);
verifyGreaterThan(testCase,reciprocalCondition,sqrt(eps));
verifyLessThan(testCase,norm(Ja*ddq+jdotdq-command),1e-10);
verifyTrue(testCase,all(isfinite([ddq;tension])));
verifyEqual(testCase,isPullOnlyFeasible,all(tension>=-sqrt(eps)));
M = plant.mass(q,plant.parameters);
h = plant.bias(q,dq,plant.parameters);
Q = plant.appliedForce(q,zeros(plant.narm,1),zeros(6,1),w,plant.parameters);
verifyLessThan(testCase,norm(M*ddq+h+Ja.'*tension-Q),1e-10);
end

function testSignedEquivalentTendon(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_signed_pair_n2"));
q = [0.01;-0.02;0.46;0.03;0.01;0.51];
tension = softarm.packActuatorInputs(plant.actuation, ...
    struct("bend_x",-2.0,"bend_y",1.5));
[tau,isFeasible] = plant.actuation.force(q,tension,plant.parameters);
Ja = plant.actuation.jacobian(q,plant.parameters);
verifyTrue(testCase,isFeasible);
verifyEqual(testCase,Ja(:,[3,6]),zeros(2,2),"AbsTol",1e-12);
verifyEqual(testCase,tau([3,6]),zeros(2,1),"AbsTol",1e-12);
verifyEqual(testCase,tau,-Ja.'*tension,"AbsTol",1e-12);
end

function testActuatorSlxArtifacts(testCase)
root = testCase.TestData.root;
libraryDirectory = fullfile(root,"matlab","simulink");
exampleDirectory = fullfile(root,"examples","simulink");
libraryFiles = [ ...
    "softarm_plant.slx", ...
    "softarm_actuator_force_block.slx", ...
    "softarm_actuator_acceleration_block.slx", ...
    "softarm_constrained_plant.slx"];
exampleFiles = [ ...
    "softarm_tendon_force_demo.slx", ...
    "softarm_tendon_acceleration_demo.slx", ...
    "softarm_flying_contact_demo.slx"];
for file = libraryFiles
verifyTrue(testCase,isfile(fullfile(libraryDirectory,file)));
end
for file = exampleFiles
verifyTrue(testCase,isfile(fullfile(exampleDirectory,file)));
end

load_system(fullfile(libraryDirectory,"softarm_plant.slx"));
cleanupPlant = onCleanup(@() close_system("softarm_plant",0));
mask = Simulink.Mask.get("softarm_plant");
verifyNotEmpty(testCase,mask);
verifyEqual(testCase,{mask.Parameters.Name},{'Bundle'});
verifyEqual(testCase,get_param("softarm_plant","ParameterArgumentNames"),'Bundle');
outputs = rootOutports("softarm_plant");
verifyEqual(testCase,outputs, ...
    ["q_out","dq_out","tip_pose","diagnostic","backbone_poses"]);

load_system(fullfile(libraryDirectory,"softarm_constrained_plant.slx"));
cleanupConstrained = onCleanup(@() close_system("softarm_constrained_plant",0));
outputs = rootOutports("softarm_constrained_plant");
verifyEqual(testCase,outputs,["q_out","dq_out","tip_pose","reaction", ...
    "is_feasible","constraint_rcond","backbone_poses"]);

load_system(fullfile(exampleDirectory,"softarm_tendon_force_demo.slx"));
cleanupForce = onCleanup(@() close_system("softarm_tendon_force_demo",0));
references = find_system("softarm_tendon_force_demo", ...
    "LookUnderMasks","all","BlockType","ModelReference");
models = string(get_param(references,"ModelName"));
verifyTrue(testCase,all(ismember( ...
    ["softarm_actuator_force_block","softarm_plant"],models)));
verifyQLog(testCase,"softarm_tendon_force_demo","Plant");
verifyEqual(testCase,get_param( ...
    "softarm_tendon_force_demo","EnablePacing"),'off');
verifyEqual(testCase,get_param( ...
    "softarm_tendon_force_demo","ReturnWorkspaceOutputs"),'off');
verifyTrue(testCase,isfile(fullfile(root,"matlab","softarm_pose_playback.m")));
verifyTrue(testCase,isfile(fullfile(root,"matlab","+softarm", ...
    "playbackBackbonePoses.m")));

load_system(fullfile(exampleDirectory,"softarm_tendon_acceleration_demo.slx"));
cleanupAcceleration = onCleanup(@() close_system( ...
    "softarm_tendon_acceleration_demo",0));
references = find_system("softarm_tendon_acceleration_demo", ...
    "LookUnderMasks","all","BlockType","ModelReference");
models = string(get_param(references,"ModelName"));
verifyTrue(testCase,all(ismember( ...
    ["softarm_actuator_acceleration_block","softarm_plant"],models)));
verifyQLog(testCase,"softarm_tendon_acceleration_demo","Plant");
verifyEqual(testCase,get_param( ...
    "softarm_tendon_acceleration_demo","ReturnWorkspaceOutputs"),'off');

load_system(fullfile(exampleDirectory,"softarm_flying_contact_demo.slx"));
cleanupContact = onCleanup(@() close_system( ...
    "softarm_flying_contact_demo",0));
verifyQLog(testCase,"softarm_flying_contact_demo","Constrained Plant");
verifyEqual(testCase,get_param( ...
    "softarm_flying_contact_demo","ReturnWorkspaceOutputs"),'off');
contactInit = string(get_param("softarm_flying_contact_demo","InitFcn"));
verifyTrue(testCase,contains(contactInit, ...
    "softarm_root=fileparts(fileparts(softarm_model_dir))"));
verifyFalse(testCase,contains(contactInit, ...
    "softarm_root=fileparts(get_param(bdroot,'FileName'))"));
set_param("softarm_flying_contact_demo","SimulationCommand","update");
verifyFalse(testCase,isfile(fullfile(root,"softarm_visualization_lib.slx")));
end

function testPoseStackWidths(testCase)
for nodeCount = 1:3
    transforms = repmat(eye(4),1,1,nodeCount);
    [translation,rpy] = softarm.poseStackToAero(transforms(:));
    verifySize(testCase,translation,[nodeCount 3]);
    verifySize(testCase,rpy,[nodeCount 3]);
end
end

function testEulerTwoSignedPairBundle(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "euler_bernoulli_ritz_two_signed_pairs_n2"));
verifyEqual(testCase,plant.nq,4);
verifyEqual(testCase,plant.actuation.count,2);
verifyEqual(testCase,string(plant.actuation.names),["pair_x","pair_y"]);
verifyTrue(testCase,isfield(plant.actuation,"acceleration"));
end

function rotation = rotationZYX(angles)
roll = angles(1); pitch = angles(2); yaw = angles(3);
cx = cos(roll); sx = sin(roll);
cy = cos(pitch); sy = sin(pitch);
cz = cos(yaw); sz = sin(yaw);
rotation = [ ...
    cz*cy, cz*sy*sx-sz*cx, cz*sy*cx+sz*sx; ...
    sz*cy, sz*sy*sx+cz*cx, sz*sy*cx-cz*sx; ...
    -sy, cy*sx, cy*cx];
end

function transform = transformZYX(position,angles)
transform = [rotationZYX(angles),reshape(position,3,1);0,0,0,1];
end

function plant = syntheticPlaybackPlant(mode,mountXYZ,mountRPY,nq,nbase)
model = struct("base_mode",mode,"segments",1, ...
    "mount_xyz",mountXYZ,"mount_rpy",mountRPY);
plant = struct("nq",nq,"nbase",nbase,"parameters",zeros(0,1), ...
    "manifest",struct("model",model),"kinematics",@identityKinematics, ...
    "centerline",@syntheticCenterline);
end

function transforms = identityKinematics(~,~)
transforms = eye(4);
end

function points = syntheticCenterline(~,~,xi)
points = zeros(3,numel(xi));
end

function q = nominalConfiguration(plant)
q = zeros(plant.nq,1);
armCoordinates = cellstr(plant.manifest.coordinates.arm);
for armIndex = 1:plant.narm
    lengthToken = regexp(armCoordinates{armIndex},'^l(\d+)$','tokens','once');
    axialToken = regexp(armCoordinates{armIndex},'^vz(\d+)$','tokens','once');
    if ~isempty(lengthToken)
        parameterName = "s"+lengthToken{1}+"_rest_length";
        parameterIndex = find(strcmp({plant.manifest.parameters.name},parameterName),1);
        q(plant.nbase+armIndex) = plant.parameters(parameterIndex);
    elseif ~isempty(axialToken)
        q(plant.nbase+armIndex) = 0;
    end
end
end

function closeIfValid(viewer)
if isvalid(viewer)
    close(viewer);
end
end

function names = rootOutports(model)
ports = find_system(model,"SearchDepth",1,"BlockType","Outport");
[~,order] = sort(str2double(get_param(ports,"Port")));
names = string(get_param(ports(order),"Name")).';
end

function verifyQLog(testCase,model,source)
logBlock = model+"/q log";
verifyEqual(testCase,get_param(logBlock,"BlockType"),'ToWorkspace');
verifyEqual(testCase,get_param(logBlock,"VariableName"),'softarm_q_log');
verifyEqual(testCase,get_param(logBlock,"SaveFormat"),'Timeseries');
verifyEqual(testCase,get_param(logBlock,"SampleTime"),'-1');
connectivity = get_param(logBlock,"PortConnectivity");
verifyEqual(testCase,string(getfullname(connectivity(1).SrcBlock)),model+"/"+source);
verifyEqual(testCase,connectivity(1).SrcPort,0);
end

function testFloatingPlaneContactBundle(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "extensible_euler_bernoulli_pcs_flying_plane_contact_n1"));
verifyEqual(testCase,plant.nbase,6);
verifyEqual(testCase,plant.narm,3);
verifyEqual(testCase,plant.constraint.family,"plane_point_contact");
q = zeros(plant.nq,1);
lengthIndex = find(strcmp({plant.manifest.parameters.name},"s1_rest_length"),1);
q(plant.nbase+3) = plant.parameters(lengthIndex);
dq = zeros(plant.nq,1);
phi = plant.constraint.value(q,plant.parameters);
verifyLessThan(testCase,abs(phi),1e-12);
[ddq,reaction,isFeasible,condition] = plant.constraint.acceleration( ...
    q,dq,zeros(plant.narm,1),zeros(6,1),zeros(6,1),zeros(1,1),plant.parameters);
A = plant.constraint.jacobian(q,plant.parameters);
gamma = plant.constraint.velocityBias(q,dq,plant.parameters);
verifyLessThan(testCase,norm(A*ddq+gamma),1e-9);
verifyTrue(testCase,all(isfinite([ddq;reaction;condition])));
verifyEqual(testCase,isFeasible,reaction>=-sqrt(eps));
M = plant.mass(q,plant.parameters);
verifyGreaterThan(testCase,norm(M(1:6,7:end),"fro"),0);
badParameters = plant.parameters;
normalIndex = find(startsWith(string({plant.manifest.parameters.name}), ...
    "constraint_plane_normal_"));
badParameters(normalIndex) = 0;
verifyError(testCase,@() plant.constraint.acceleration( ...
    q,dq,zeros(plant.narm,1),zeros(6,1),zeros(6,1),zeros(1,1),badParameters), ...
    "softarm:InvalidPlaneNormal");
end
