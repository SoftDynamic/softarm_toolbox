function build_softarm_actuator_models
%BUILD_SOFTARM_ACTUATOR_MODELS Build reusable actuator models and GUI demos.
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
libraryDirectory = fullfile(repositoryRoot,"matlab","simulink");
exampleDirectory = fullfile(repositoryRoot,"examples","simulink");
if ~isfolder(libraryDirectory), mkdir(libraryDirectory); end
if ~isfolder(exampleDirectory), mkdir(exampleDirectory); end
if any(strcmp(strsplit(path,pathsep),libraryDirectory)), rmpath(libraryDirectory); end
buildForceActuator(fullfile(libraryDirectory,"softarm_actuator_force_block.slx"));
buildAccelerationActuator(fullfile(libraryDirectory,"softarm_actuator_acceleration_block.slx"));
addpath(libraryDirectory);
buildActuationDemo(fullfile(exampleDirectory,"softarm_tendon_actuation_demo.slx"));
end

function code = actuatorInit(requireAcceleration)
condition = "assert(~isempty(softarm_loaded.actuation),'softarm:MissingActuation','Selected bundle has no actuator mapping.');";
if requireAcceleration
    condition = condition + "assert(isfield(softarm_loaded.actuation,'acceleration'),'softarm:MissingAcceleration','Selected bundle has no strict acceleration mapping.');";
end
code = "softarm_model_dir=fileparts(get_param(bdroot,'FileName'));" + ...
    "softarm_root=fileparts(fileparts(softarm_model_dir));" + ...
    "addpath(fullfile(softarm_root,'matlab'));" + ...
    "addpath(fullfile(softarm_root,'matlab','simulink'));" + ...
    "if evalin('base','exist(''softarm_bundle'',''var'')'),softarm_selected_bundle=evalin('base','softarm_bundle');else,softarm_selected_bundle=fullfile(softarm_root,'examples','generated','extensible_euler_bernoulli_pcs_three_tendon_n2');end;" + ...
    "softarm_loaded=softarm.initModel(string(softarm_selected_bundle));" + condition;
end

function code = actuationDemoInit
code = "softarm_model_dir=fileparts(get_param(bdroot,'FileName'));" + ...
    "softarm_root=fileparts(fileparts(softarm_model_dir));" + ...
    "addpath(fullfile(softarm_root,'matlab'));" + ...
    "addpath(fullfile(softarm_root,'matlab','simulink'));" + ...
    "softarm_plant_block=string(bdroot)+""/Plant"";" + ...
    "softarm_actuator_block=string(bdroot)+""/Actuator mode"";" + ...
    "softarm_mode=string(get_param(softarm_actuator_block,'Mode'));" + ...
    "set_param(softarm_actuator_block,'LabelModeActiveChoice',char(softarm_mode));" + ...
    "if softarm_mode==""Acceleration"",softarm_requirement=""strict"";else,softarm_requirement=""force"";end;" + ...
    "softarm.applyPlantMask(softarm_plant_block,softarm_requirement);";
end

function prepare(model, initCode, stopTime)
if bdIsLoaded(model), close_system(model,0); end
new_system(model);
set_param(model,"Solver","ode15s","StopTime",stopTime,"InitFcn",initCode);
set_param(model,"PreLoadFcn", ...
    "softarm_model_dir=fileparts(get_param(bdroot,'FileName'));" + ...
    "softarm_root=fileparts(fileparts(softarm_model_dir));" + ...
    "addpath(fullfile(softarm_root,'matlab'));" + ...
    "addpath(fullfile(softarm_root,'matlab','simulink'));");
end

function addIn(model, name, position, dimension)
add_block("simulink/Sources/In1",model+"/"+name, ...
    "Position",position,"PortDimensions",dimension);
end

function addOut(model, name, position)
add_block("simulink/Sinks/Out1",model+"/"+name,"Position",position);
end

function addQLog(model,source,position)
add_block("simulink/Sinks/To Workspace",model+"/q log", ...
    "Position",position,"VariableName","softarm_q_log", ...
    "SaveFormat","Timeseries","SampleTime","-1", ...
    "Decimation","1","MaxDataPoints","inf");
add_line(model,source,"q log/1","autorouting","on");
end

function buildForceActuator(outputPath)
model = "softarm_actuator_force_block";
prepare(model,actuatorInit(false),"1");
addIn(model,"q",[30 65 60 85],"softarm_nq");
addIn(model,"action",[30 125 60 145],"softarm_nu");
add_block("simulink/Sources/Constant",model+"/parameters", ...
    "Position",[30 185 90 215],"Value","softarm_p");
add_block("simulink/User-Defined Functions/MATLAB Function",model+"/Tendon force map", ...
    "Position",[155 65 310 205]);
chart = find(sfroot,"-isa","Stateflow.EMChart","Path",model+"/Tendon force map");
chart.Script = sprintf("function [tau,isFeasible] = fcn(q,action,p)\n%%#codegen\n[tau,isFeasible]=softarm_actuator_force(q,action,p);\nend\n");
addOut(model,"tau",[390 90 420 110]);
addOut(model,"is_feasible",[390 155 420 175]);
add_line(model,"q/1","Tendon force map/1");
add_line(model,"action/1","Tendon force map/2");
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
addIn(model,"action",[30 130 60 150],"softarm_nu");
addIn(model,"tau_arm_external",[30 175 60 195],"softarm_narm");
addIn(model,"vehicle_wrench",[30 220 60 240],"6");
addIn(model,"tip_wrench",[30 265 60 285],"6");
add_block("simulink/Sources/Constant",model+"/parameters", ...
    "Position",[30 310 90 340],"Value","softarm_p");
add_block("simulink/User-Defined Functions/MATLAB Function",model+"/Acceleration constraint", ...
    "Position",[155 45 365 335]);
chart = find(sfroot,"-isa","Stateflow.EMChart","Path",model+"/Acceleration constraint");
chart.Script = sprintf("function [tau,tension,isFeasible,condition,ddq] = fcn(q,dq,action,tauArmExternal,wVehicle,wTip,p)\n%%#codegen\n[ddq,tension,isFeasible,condition]=softarm_actuator_acceleration(q,dq,action,tauArmExternal,wVehicle,wTip,p);\nJa=softarm_actuator_jacobian(q,p);\ntau=tauArmExternal-Ja.'*tension;\nend\n");
outputs = ["tau","tension","is_feasible","constraint_rcond","constrained_ddq"];
for index = 1:numel(outputs)
    addOut(model,outputs(index),[440 35+50*index 470 55+50*index]);
end
for index = 1:6
    add_line(model,modelPortName(index)+"/1","Acceleration constraint/"+index);
end
add_line(model,"parameters/1","Acceleration constraint/7");
for index = 1:5
    add_line(model,"Acceleration constraint/"+index,outputs(index)+"/1");
end
save_system(model,outputPath);
close_system(model,0);
end

function name = modelPortName(index)
names = ["q","dq","action","tau_arm_external","vehicle_wrench","tip_wrench"];
name = names(index);
end

function addReference(model, name, referencedModel, position)
add_block("built-in/ModelReference",model+"/"+name, ...
    "ModelName",referencedModel,"Position",position);
end

function buildActuationDemo(outputPath)
model = "softarm_tendon_actuation_demo";
prepare(model,actuationDemoInit,"0.5");
save_system(model,outputPath);
set_param(model,"EnablePacing","off","ReturnWorkspaceOutputs","off");
add_block("simulink/Sources/Constant",model+"/Action command", ...
    "Position",[20 60 135 90],"Value","softarm_u0");
add_block("simulink/Sources/Constant",model+"/Additional tau", ...
    "Position",[20 115 135 145],"Value","zeros(softarm_narm,1)");
add_block("simulink/Sources/Constant",model+"/Vehicle wrench", ...
    "Position",[20 170 135 200],"Value","zeros(6,1)");
add_block("simulink/Sources/Constant",model+"/Tip wrench", ...
    "Position",[20 225 135 255],"Value","zeros(6,1)");
addActuatorVariant(model,"Actuator mode",[205 35 435 275]);
addReference(model,"Plant","softarm_plant",[540 55 695 200]);
set_param(model+"/Plant","Bundle","'../generated/extensible_euler_bernoulli_pcs_three_tendon_n2'");
add_block("simulink/Discrete/Memory",model+"/q feedback memory", ...
    "Position",[540 250 575 280],"InitialCondition","softarm_x0(1:softarm_nq)");
add_block("simulink/Discrete/Memory",model+"/dq feedback memory", ...
    "Position",[605 300 640 330],"InitialCondition","zeros(softarm_nq,1)");
addOut(model,"q",[790 45 820 65]);
addOut(model,"dq",[790 85 820 105]);
addOut(model,"tip_pose",[790 125 820 145]);
addOut(model,"mass_rcond",[790 165 820 185]);
addOut(model,"tension",[470 75 500 95]);
addOut(model,"tension_feasible",[470 120 500 140]);
addOut(model,"constraint_rcond",[470 165 500 185]);
add_block("simulink/Sinks/Scope",model+"/State scope", ...
    "Position",[875 60 910 120],"NumInputPorts","2");
set_param(model,"SimulationCommand","update");
addQLog(model,"Plant/1",[875 145 975 175]);
add_line(model,"Plant/1","q feedback memory/1","autorouting","on");
add_line(model,"q feedback memory/1","Actuator mode/1","autorouting","on");
add_line(model,"Plant/2","dq feedback memory/1","autorouting","on");
add_line(model,"dq feedback memory/1","Actuator mode/2","autorouting","on");
add_line(model,"Action command/1","Actuator mode/3");
add_line(model,"Additional tau/1","Actuator mode/4");
add_line(model,"Vehicle wrench/1","Actuator mode/5");
add_line(model,"Tip wrench/1","Actuator mode/6");
add_line(model,"Actuator mode/1","Plant/1");
add_line(model,"Vehicle wrench/1","Plant/2","autorouting","on");
add_line(model,"Tip wrench/1","Plant/3","autorouting","on");
add_line(model,"Actuator mode/2","tension/1");
add_line(model,"Actuator mode/3","tension_feasible/1");
add_line(model,"Actuator mode/4","constraint_rcond/1");
destinations = ["q","dq","tip_pose","mass_rcond"];
for index = 1:4
    add_line(model,"Plant/"+index,destinations(index)+"/1");
end
add_line(model,"Plant/1","State scope/1");
add_line(model,"Plant/2","State scope/2");
save_system(model,outputPath);
close_system(model,0);
end

function addActuatorVariant(model, name, position)
variant = model+"/"+name;
source = sprintf("simulink/Ports &\nSubsystems/Variant Subsystem");
add_block(source,variant,"Position",position);
delete_block(variant+"/In1");
delete_block(variant+"/Out1");
addChoicePorts(variant);
set_param(variant+"/Subsystem","Name","Force");
set_param(variant+"/Subsystem1","Name","Acceleration");
set_param(variant+"/Force","VariantControl","Force");
set_param(variant+"/Acceleration","VariantControl","Acceleration");
set_param(variant,"VariantControlMode","label", ...
    "LabelModeActiveChoice","Force");
buildForceChoice(variant+"/Force");
buildAccelerationChoice(variant+"/Acceleration");

mask = Simulink.Mask.create(variant);
mask.Description = "Select tendon tension or strict tendon-coordinate acceleration actuation.";
modeParameter = mask.addParameter("Type","popup","Name","Mode", ...
    "Prompt","Actuation mode","TypeOptions",{'Force','Acceleration'}, ...
    "Value","Force");
modeParameter.Callback = ...
    "softarm_mode=string(get_param(gcb,'Mode'));" + ...
    "set_param(gcb,'LabelModeActiveChoice',char(softarm_mode));" + ...
    "if softarm_mode==""Acceleration"",softarm_requirement=""strict"";else,softarm_requirement=""force"";end;" + ...
    "softarm.applyPlantMask(string(bdroot(gcb))+""/Plant"",softarm_requirement);";
end

function buildForceChoice(choice)
Simulink.SubSystem.deleteContents(choice);
addChoicePorts(choice);
addReference(choice,"Force actuator","softarm_actuator_force_block",[175 35 345 125]);
add_block("simulink/Sources/Constant",choice+"/constraint unavailable", ...
    "Position",[360 205 410 235],"Value","NaN");
add_line(choice,"q/1","Force actuator/1");
add_line(choice,"action/1","Force actuator/2");
add_line(choice,"Force actuator/1","tau/1");
add_line(choice,"action/1","tension/1","autorouting","on");
add_line(choice,"Force actuator/2","is_feasible/1");
add_line(choice,"constraint unavailable/1","constraint_rcond/1");
addTerminator(choice,"unused dq","dq/1",[105 55 125 75]);
addTerminator(choice,"unused tau","tau_arm_external/1",[105 135 125 155]);
addTerminator(choice,"unused vehicle wrench","vehicle_wrench/1",[105 175 125 195]);
addTerminator(choice,"unused tip wrench","tip_wrench/1",[105 215 125 235]);
end

function buildAccelerationChoice(choice)
Simulink.SubSystem.deleteContents(choice);
addChoicePorts(choice);
addReference(choice,"Acceleration actuator", ...
    "softarm_actuator_acceleration_block",[175 25 365 245]);
for index = 1:6
    add_line(choice,modelPortName(index)+"/1","Acceleration actuator/"+index);
end
outputs = ["tau","tension","is_feasible","constraint_rcond"];
for index = 1:numel(outputs)
    add_line(choice,"Acceleration actuator/"+index,outputs(index)+"/1");
end
addTerminator(choice,"unused constrained ddq","Acceleration actuator/5", ...
    [395 225 415 245]);
end

function addChoicePorts(choice)
inputs = modelPortName(1:6);
dimensions = ["softarm_nq","softarm_nq","softarm_nu", ...
    "softarm_narm","6","6"];
for index = 1:numel(inputs)
    addIn(choice,inputs(index),[25 20+40*index 55 40+40*index], ...
        dimensions(index));
end
outputs = ["tau","tension","is_feasible","constraint_rcond"];
for index = 1:numel(outputs)
    addOut(choice,outputs(index),[455 45+50*index 485 65+50*index]);
end
end

function addTerminator(model, name, source, position)
add_block("simulink/Sinks/Terminator",model+"/"+name,"Position",position);
add_line(model,source,name+"/1","autorouting","on");
end
