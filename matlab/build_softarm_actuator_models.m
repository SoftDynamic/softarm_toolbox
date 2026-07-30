function build_softarm_actuator_models(outputDirectory)
%BUILD_SOFTARM_ACTUATOR_MODELS Build reusable actuator models and GUI demos.
arguments
    outputDirectory (1,1) string = string(fileparts(fileparts(mfilename("fullpath"))))
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
buildForceActuator(fullfile(outputDirectory,"softarm_actuator_force_block.slx"));
buildAccelerationActuator(fullfile(outputDirectory,"softarm_actuator_acceleration_block.slx"));
buildForceDemo(fullfile(outputDirectory,"softarm_tendon_force_demo.slx"));
buildAccelerationDemo(fullfile(outputDirectory,"softarm_tendon_acceleration_demo.slx"));
end

function code = actuatorInit(requireAcceleration)
condition = "assert(~isempty(softarm_loaded.actuation),'softarm:MissingActuation','Selected bundle has no actuator mapping.');";
if requireAcceleration
    condition = condition + "assert(isfield(softarm_loaded.actuation,'acceleration'),'softarm:MissingAcceleration','Selected bundle has no strict acceleration mapping.');";
end
code = "softarm_root=fileparts(get_param(bdroot,'FileName'));addpath(fullfile(softarm_root,'matlab'));" + ...
    "if evalin('base','exist(''softarm_bundle'',''var'')'),softarm_selected_bundle=evalin('base','softarm_bundle');else,softarm_selected_bundle=fullfile(softarm_root,'examples','generated','pcc_three_tendon_extensible_n2');end;" + ...
    "softarm_loaded=softarm.initModel(string(softarm_selected_bundle));" + condition;
end

function code = demoInit(requirement)
code = "softarm_root=fileparts(get_param(bdroot,'FileName'));addpath(fullfile(softarm_root,'matlab'));" + ...
    "softarm_plant_block=string(bdroot)+""/Plant"";" + ...
    "softarm.applyPlantMask(softarm_plant_block,'" + requirement + "');";
end

function prepare(model, initCode, stopTime)
if bdIsLoaded(model), close_system(model,0); end
new_system(model);
set_param(model,"Solver","ode15s","StopTime",stopTime,"InitFcn",initCode);
set_param(model,"PreLoadFcn","softarm_root=fileparts(get_param(bdroot,'FileName'));addpath(softarm_root);addpath(fullfile(softarm_root,'matlab'));");
end

function addIn(model, name, position, dimension)
add_block("simulink/Sources/In1",model+"/"+name, ...
    "Position",position,"PortDimensions",dimension);
end

function addOut(model, name, position)
add_block("simulink/Sinks/Out1",model+"/"+name,"Position",position);
end

function buildForceActuator(outputPath)
model = "softarm_actuator_force_block";
prepare(model,actuatorInit(false),"1");
addIn(model,"q",[30 65 60 85],"softarm_nq");
addIn(model,"tension",[30 125 60 145],"softarm_nu");
add_block("simulink/Sources/Constant",model+"/parameters", ...
    "Position",[30 185 90 215],"Value","softarm_p");
add_block("simulink/User-Defined Functions/MATLAB Function",model+"/Tendon force map", ...
    "Position",[155 65 310 205]);
chart = find(sfroot,"-isa","Stateflow.EMChart","Path",model+"/Tendon force map");
chart.Script = sprintf("function [tau,isFeasible] = fcn(q,tension,p)\n%%#codegen\n[tau,isFeasible]=softarm_actuator_force(q,tension,p);\nend\n");
addOut(model,"tau",[390 90 420 110]);
addOut(model,"is_feasible",[390 155 420 175]);
add_line(model,"q/1","Tendon force map/1");
add_line(model,"tension/1","Tendon force map/2");
add_line(model,"parameters/1","Tendon force map/3");
add_line(model,"Tendon force map/1","tau/1");
add_line(model,"Tendon force map/2","is_feasible/1");
save_system(model,outputPath);
close_system(model,0);
end

function buildAccelerationActuator(outputPath)
model = "softarm_actuator_acceleration_block";
prepare(model,actuatorInit(true),"1");
addIn(model,"q",[30 40 60 60],"softarm_nq");
addIn(model,"dq",[30 85 60 105],"softarm_nq");
addIn(model,"coordinate_acceleration",[30 130 60 150],"softarm_nu");
addIn(model,"tau_external",[30 175 60 195],"softarm_nq");
addIn(model,"w",[30 220 60 240],"6");
add_block("simulink/Sources/Constant",model+"/parameters", ...
    "Position",[30 265 90 295],"Value","softarm_p");
add_block("simulink/User-Defined Functions/MATLAB Function",model+"/Acceleration constraint", ...
    "Position",[155 45 350 290]);
chart = find(sfroot,"-isa","Stateflow.EMChart","Path",model+"/Acceleration constraint");
chart.Script = sprintf("function [tau,tension,isFeasible,condition,ddq] = fcn(q,dq,ydd,tauExternal,w,p)\n%%#codegen\n[ddq,tension,isFeasible,condition]=softarm_actuator_acceleration(q,dq,ydd,tauExternal,w,p);\nJa=softarm_actuator_jacobian(q,p);\ntau=tauExternal-Ja.'*tension;\nend\n");
outputs = ["tau","tension","is_feasible","constraint_rcond","constrained_ddq"];
for index = 1:numel(outputs)
    addOut(model,outputs(index),[440 35+50*index 470 55+50*index]);
end
for index = 1:5
    add_line(model,modelPortName(index)+"/1","Acceleration constraint/"+index);
end
add_line(model,"parameters/1","Acceleration constraint/6");
for index = 1:5
    add_line(model,"Acceleration constraint/"+index,outputs(index)+"/1");
end
save_system(model,outputPath);
close_system(model,0);
end

function name = modelPortName(index)
names = ["q","dq","coordinate_acceleration","tau_external","w"];
name = names(index);
end

function addReference(model, name, referencedModel, position)
add_block("built-in/ModelReference",model+"/"+name, ...
    "ModelName",referencedModel,"Position",position);
end

function buildForceDemo(outputPath)
model = "softarm_tendon_force_demo";
prepare(model,demoInit("force"),"0.5");
add_block("simulink/Sources/Constant",model+"/Tension command", ...
    "Position",[30 65 120 95],"Value","softarm_u0");
add_block("simulink/Sources/Constant",model+"/External wrench", ...
    "Position",[360 185 450 215],"Value","zeros(6,1)");
addReference(model,"Actuator","softarm_actuator_force_block",[190 40 350 120]);
addReference(model,"Plant","softarm_plant",[510 65 665 210]);
set_param(model+"/Plant","Bundle","'examples/generated/pcc_three_tendon_extensible_n2'");
add_block("simulink/Discrete/Memory",model+"/q feedback memory", ...
    "Position",[510 260 545 290],"InitialCondition","softarm_x0(1:softarm_nq)");
addOut(model,"q",[760 55 790 75]);
addOut(model,"dq",[760 100 790 120]);
addOut(model,"tip_pose",[760 145 790 165]);
addOut(model,"mass_rcond",[760 190 790 210]);
addOut(model,"tension_feasible",[410 105 440 125]);
add_block("simulink/Sinks/Scope",model+"/State scope", ...
    "Position",[850 70 885 130],"NumInputPorts","2");
add_line(model,"Tension command/1","Actuator/2");
add_line(model,"Actuator/1","Plant/1");
add_line(model,"Actuator/2","tension_feasible/1");
add_line(model,"External wrench/1","Plant/2");
add_line(model,"Plant/1","q feedback memory/1","autorouting","on");
add_line(model,"q feedback memory/1","Actuator/1","autorouting","on");
for index = 1:4
    destinations = ["q","dq","tip_pose","mass_rcond"];
    add_line(model,"Plant/"+index,destinations(index)+"/1");
end
add_line(model,"Plant/1","State scope/1");
add_line(model,"Plant/2","State scope/2");
save_system(model,outputPath);
close_system(model,0);
end

function buildAccelerationDemo(outputPath)
model = "softarm_tendon_acceleration_demo";
prepare(model,demoInit("strict"),"0.5");
add_block("simulink/Sources/Constant",model+"/Tendon acceleration", ...
    "Position",[20 60 135 90],"Value","zeros(softarm_nu,1)");
add_block("simulink/Sources/Constant",model+"/Additional tau", ...
    "Position",[20 115 135 145],"Value","zeros(softarm_nq,1)");
add_block("simulink/Sources/Constant",model+"/External wrench", ...
    "Position",[20 170 135 200],"Value","zeros(6,1)");
addReference(model,"Acceleration actuator","softarm_actuator_acceleration_block",[210 35 405 205]);
addReference(model,"Plant","softarm_plant",[520 55 675 200]);
set_param(model+"/Plant","Bundle","'examples/generated/pcc_three_tendon_extensible_n2'");
add_block("simulink/Discrete/Memory",model+"/q feedback memory", ...
    "Position",[520 245 555 275],"InitialCondition","softarm_x0(1:softarm_nq)");
add_block("simulink/Discrete/Memory",model+"/dq feedback memory", ...
    "Position",[585 285 620 315],"InitialCondition","zeros(softarm_nq,1)");
addOut(model,"q",[770 45 800 65]);
addOut(model,"dq",[770 85 800 105]);
addOut(model,"tip_pose",[770 125 800 145]);
addOut(model,"mass_rcond",[770 165 800 185]);
addOut(model,"tension",[455 90 485 110]);
addOut(model,"tension_feasible",[455 125 485 145]);
addOut(model,"constraint_rcond",[455 160 485 180]);
add_block("simulink/Sinks/Scope",model+"/State scope", ...
    "Position",[855 60 890 120],"NumInputPorts","2");
add_line(model,"Plant/1","q feedback memory/1","autorouting","on");
add_line(model,"q feedback memory/1","Acceleration actuator/1","autorouting","on");
add_line(model,"Plant/2","dq feedback memory/1","autorouting","on");
add_line(model,"dq feedback memory/1","Acceleration actuator/2","autorouting","on");
add_line(model,"Tendon acceleration/1","Acceleration actuator/3");
add_line(model,"Additional tau/1","Acceleration actuator/4");
add_line(model,"External wrench/1","Acceleration actuator/5");
add_line(model,"Acceleration actuator/1","Plant/1");
add_line(model,"External wrench/1","Plant/2","autorouting","on");
add_line(model,"Acceleration actuator/2","tension/1");
add_line(model,"Acceleration actuator/3","tension_feasible/1");
add_line(model,"Acceleration actuator/4","constraint_rcond/1");
destinations = ["q","dq","tip_pose","mass_rcond"];
for index = 1:4
    add_line(model,"Plant/"+index,destinations(index)+"/1");
end
add_line(model,"Plant/1","State scope/1");
add_line(model,"Plant/2","State scope/2");
save_system(model,outputPath);
close_system(model,0);
end
