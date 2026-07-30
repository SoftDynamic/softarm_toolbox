function tests = test_generated
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
testCase.TestData.root = root;
end

function testReferenceBundles(testCase)
names = ["pcc_lumped_n2","pcc_distributed_n2","euler_ritz_n2"];
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

function testLinearization(testCase)
bundle = fullfile(testCase.TestData.root,"examples","generated","euler_ritz_n2");
plant = softarm.loadModel(bundle);
x = zeros(2*plant.nq,1); tau = zeros(plant.narm,1);
[A,Barm,Bvehicle,Btip] = softarm.linearize(plant,x,tau);
verifySize(testCase,A,[2*plant.nq,2*plant.nq]);
verifySize(testCase,Barm,[2*plant.nq,plant.narm]);
verifySize(testCase,Bvehicle,[2*plant.nq,6]);
verifySize(testCase,Btip,[2*plant.nq,6]);
verifyTrue(testCase,all(isfinite([A(:);Barm(:);Bvehicle(:);Btip(:)])));
end

function testThreeTendonForceAndAxialMode(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "pcc_three_tendon_extensible_n2"));
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
    "pcc_three_tendon_extensible_n2"));
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
plant = softarm.loadModel(fullfile(root,"examples","generated","pcc_signed_pair_n2"));
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
files = [ ...
    "softarm_plant.slx", ...
    "softarm_actuator_force_block.slx", ...
    "softarm_actuator_acceleration_block.slx", ...
    "softarm_tendon_force_demo.slx", ...
    "softarm_tendon_acceleration_demo.slx", ...
    "softarm_constrained_plant.slx", ...
    "softarm_flying_contact_demo.slx"];
for file = files
    verifyTrue(testCase,isfile(fullfile(root,file)));
end

load_system(fullfile(root,"softarm_plant.slx"));
cleanupPlant = onCleanup(@() close_system("softarm_plant",0)); %#ok<NASGU>
mask = Simulink.Mask.get("softarm_plant");
verifyNotEmpty(testCase,mask);
verifyEqual(testCase,{mask.Parameters.Name},{'Bundle'});
verifyEqual(testCase,get_param("softarm_plant","ParameterArgumentNames"),'Bundle');

load_system(fullfile(root,"softarm_tendon_force_demo.slx"));
cleanupForce = onCleanup(@() close_system("softarm_tendon_force_demo",0)); %#ok<NASGU>
references = find_system("softarm_tendon_force_demo", ...
    "LookUnderMasks","all","BlockType","ModelReference");
models = string(get_param(references,"ModelName"));
verifyTrue(testCase,all(ismember( ...
    ["softarm_actuator_force_block","softarm_plant"],models)));

load_system(fullfile(root,"softarm_tendon_acceleration_demo.slx"));
cleanupAcceleration = onCleanup(@() close_system( ...
    "softarm_tendon_acceleration_demo",0)); %#ok<NASGU>
references = find_system("softarm_tendon_acceleration_demo", ...
    "LookUnderMasks","all","BlockType","ModelReference");
models = string(get_param(references,"ModelName"));
verifyTrue(testCase,all(ismember( ...
    ["softarm_actuator_acceleration_block","softarm_plant"],models)));
end

function testEulerTwoSignedPairBundle(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "euler_two_signed_pairs_n2"));
verifyEqual(testCase,plant.nq,4);
verifyEqual(testCase,plant.actuation.count,2);
verifyEqual(testCase,string(plant.actuation.names),["pair_x","pair_y"]);
verifyTrue(testCase,isfield(plant.actuation,"acceleration"));
end

function testFloatingPlaneContactBundle(testCase)
root = testCase.TestData.root;
plant = softarm.loadModel(fullfile(root,"examples","generated", ...
    "pcc_flying_plane_contact_n1"));
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
