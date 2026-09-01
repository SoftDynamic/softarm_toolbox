function result = verify_recursive_parity(symbolicBundle,recursiveBundle,q,dq,ddq)
%VERIFY_RECURSIVE_PARITY Compare generated symbolic and recursive bundles.
arguments
    symbolicBundle (1,1) string
    recursiveBundle (1,1) string
    q (:,1) double
    dq (:,1) double
    ddq (:,1) double
end
symbolic = evaluateBundle(symbolicBundle,q,dq,ddq);
recursive = evaluateBundle(recursiveBundle,q,dq,ddq);
result.massError = norm(recursive.mass-symbolic.mass,"fro");
result.biasError = norm(recursive.bias-symbolic.bias,inf);
result.inverseDynamicsError = norm( ...
    recursive.inverseDynamics-symbolic.inverseDynamics,inf);
result.kinematicsError = norm( ...
    recursive.kinematics(:)-symbolic.kinematics(:),inf);
result.constraintValueError = compareOptional( ...
    recursive,"constraintValue",symbolic,"constraintValue");
result.constraintJacobianError = compareOptional( ...
    recursive,"constraintJacobian",symbolic,"constraintJacobian");
result.constraintBiasError = compareOptional( ...
    recursive,"constraintBias",symbolic,"constraintBias");
result.constraintReactionMapError = compareOptional( ...
    recursive,"constraintReactionMap",symbolic,"constraintReactionMap");
assert(result.massError<1e-8,"softarm:RecursiveMassMismatch", ...
    "Recursive and symbolic mass matrices differ.");
assert(result.biasError<1e-8,"softarm:RecursiveBiasMismatch", ...
    "Recursive and symbolic bias vectors differ.");
assert(result.inverseDynamicsError<1e-8,"softarm:RecursiveIDMismatch", ...
    "Recursive and symbolic inverse dynamics differ.");
assert(result.kinematicsError<1e-9,"softarm:RecursiveKinematicsMismatch", ...
    "Recursive and symbolic kinematics differ.");
end

function values = evaluateBundle(bundle,q,dq,ddq)
manifest = jsondecode(fileread(fullfile(bundle,"manifest.json")));
p = [manifest.parameters.default].';
addpath(bundle);
cleanup = onCleanup(@() removeBundle(bundle)); %#ok<NASGU>
values.mass = softarm_mass(q,p);
values.bias = softarm_bias(q,dq,p);
if isfile(fullfile(bundle,"softarm_inverse_dynamics.m"))
    values.inverseDynamics = softarm_inverse_dynamics(q,dq,ddq,p);
else
    values.inverseDynamics = values.mass*ddq+values.bias;
end
values.kinematics = softarm_kinematics(q,p);
if ~isempty(manifest.constraint)
    values.constraintValue = softarm_constraint_value(q,p);
    values.constraintJacobian = softarm_constraint_jacobian(q,p);
    values.constraintBias = softarm_constraint_velocity_bias(q,dq,p);
    values.constraintReactionMap = softarm_constraint_reaction_map(q,dq,p);
end
end

function removeBundle(bundle)
rmpath(bundle);
clear softarm_mass softarm_bias softarm_inverse_dynamics softarm_kinematics
clear softarm_constraint_value softarm_constraint_jacobian
clear softarm_constraint_velocity_bias softarm_constraint_reaction_map
clear softarm_constraint_stabilization
end

function error = compareOptional(left,leftField,right,rightField)
if isfield(left,leftField) && isfield(right,rightField)
    error = norm(left.(leftField)-right.(rightField),inf);
else
    error = NaN;
end
end
