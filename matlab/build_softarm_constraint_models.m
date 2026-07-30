function build_softarm_constraint_models(outputDirectory)
%BUILD_SOFTARM_CONSTRAINT_MODELS Build the generic constrained Plant and demo.
arguments
    outputDirectory (1,1) string = string(fileparts(fileparts(mfilename("fullpath"))))
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
buildConstrainedPlant(fullfile(outputDirectory,"softarm_constrained_plant.slx"));
buildContactDemo(fullfile(outputDirectory,"softarm_flying_contact_demo.slx"));
end

function prepare(model, stopTime)
if bdIsLoaded(model), close_system(model,0); end
new_system(model);
set_param(model,"Solver","ode15s","StopTime",stopTime);
set_param(model,"PreLoadFcn", ...
    "softarm_root=fileparts(get_param(bdroot,'FileName'));addpath(softarm_root);addpath(fullfile(softarm_root,'matlab'));");
end

function addIn(model,name,position,dimension)
add_block("simulink/Sources/In1",model+"/"+name,"Position",position, ...
    "PortDimensions",dimension);
end

function addOut(model,name,position)
add_block("simulink/Sinks/Out1",model+"/"+name,"Position",position);
end

function buildConstrainedPlant(outputPath)
model = "softarm_constrained_plant";
prepare(model,"5");
defaultBundle = "examples/generated/pcc_flying_plane_contact_n1";
modelWorkspace = get_param(model,"ModelWorkspace");
assignin(modelWorkspace,"Bundle",char(defaultBundle));
set_param(model,"ParameterArgumentNames","Bundle");
mask = Simulink.Mask.create(model);
mask.Description = "Generated floating/fixed SoftArm Plant with an active acceleration constraint.";
parameter = mask.getParameter("Bundle");
parameter.Prompt = "Bundle directory";
parameter.Callback = "if string(get_param(gcb,'Type'))~=""block_diagram"",softarm.applyPlantMask(string(gcb),'constraint');end;";

addIn(model,"tau_arm",[25 45 55 65],"softarm_narm");
addIn(model,"vehicle_wrench",[25 90 55 110],"6");
addIn(model,"tip_wrench",[25 135 55 155],"6");
addIn(model,"constraint_acceleration",[25 180 55 200],"softarm_nc");
add_block("simulink/Sources/Constant",model+"/parameters", ...
    "Position",[25 230 85 260],"Value","softarm_p");
add_block("simulink/Continuous/Integrator",model+"/state", ...
    "Position",[365 55 395 85],"InitialCondition","softarm_x0");
add_block("simulink/Signal Routing/Demux",model+"/State split", ...
    "Position",[450 40 455 130],"Outputs","[softarm_nq softarm_nq]");

add_block("simulink/User-Defined Functions/MATLAB Function",model+"/Constrained RHS", ...
    "Position",[145 35 315 270]);
chart = find(sfroot,"-isa","Stateflow.EMChart","Path",model+"/Constrained RHS");
chart.Script = sprintf( ...
    "function [dx,reaction,isFeasible,condition] = fcn(x,tauArm,wVehicle,wTip,ac,p)\n" + ...
    "%%#codegen\nn=numel(x)/2; q=x(1:n); dq=x(n+1:end);\n" + ...
    "[ddq,reaction,isFeasible,condition]=softarm_constraint_acceleration(q,dq,tauArm,wVehicle,wTip,ac,p);\n" + ...
    "dx=[dq;ddq];\nend\n");

add_block("simulink/User-Defined Functions/MATLAB Function",model+"/Tip pose", ...
    "Position",[535 155 645 215]);
chart = find(sfroot,"-isa","Stateflow.EMChart","Path",model+"/Tip pose");
chart.Script = sprintf("function pose = fcn(q,p)\n%%#codegen\nH=softarm_kinematics(q,p); pose=reshape(H(:,:,end),16,1);\nend\n");

addOut(model,"q_out",[720 35 750 55]);
addOut(model,"dq_out",[720 80 750 100]);
addOut(model,"tip_pose",[720 165 750 185]);
addOut(model,"reaction",[370 165 400 185]);
addOut(model,"is_feasible",[370 210 400 230]);
addOut(model,"constraint_rcond",[370 255 400 275]);

add_line(model,"state/1","State split/1");
add_line(model,"state/1","Constrained RHS/1");
add_line(model,"tau_arm/1","Constrained RHS/2");
add_line(model,"vehicle_wrench/1","Constrained RHS/3");
add_line(model,"tip_wrench/1","Constrained RHS/4");
add_line(model,"constraint_acceleration/1","Constrained RHS/5");
add_line(model,"parameters/1","Constrained RHS/6");
add_line(model,"Constrained RHS/1","state/1");
add_line(model,"Constrained RHS/2","reaction/1");
add_line(model,"Constrained RHS/3","is_feasible/1");
add_line(model,"Constrained RHS/4","constraint_rcond/1");
add_line(model,"State split/1","q_out/1");
add_line(model,"State split/2","dq_out/1");
add_line(model,"State split/1","Tip pose/1");
add_line(model,"parameters/1","Tip pose/2");
add_line(model,"Tip pose/1","tip_pose/1");
save_system(model,outputPath);
close_system(model,0);
end

function buildContactDemo(outputPath)
model = "softarm_flying_contact_demo";
prepare(model,"1");
set_param(model,"InitFcn", ...
    "softarm_root=fileparts(get_param(bdroot,'FileName'));addpath(fullfile(softarm_root,'matlab'));softarm.initModel(fullfile(softarm_root,'examples','generated','pcc_flying_plane_contact_n1'));");
add_block("simulink/Sources/Constant",model+"/Arm generalized force", ...
    "Position",[25 45 135 75],"Value","zeros(softarm_narm,1)");
add_block("simulink/Sources/Constant",model+"/Vehicle wrench", ...
    "Position",[25 95 135 125],"Value","zeros(6,1)");
add_block("simulink/Sources/Constant",model+"/Tip wrench", ...
    "Position",[25 145 135 175],"Value","zeros(6,1)");
add_block("simulink/Sources/Constant",model+"/Constraint acceleration", ...
    "Position",[25 195 135 225],"Value","zeros(softarm_nc,1)");
add_block("built-in/ModelReference",model+"/Constrained Plant", ...
    "ModelName","softarm_constrained_plant","Position",[230 45 420 230]);
set_param(model+"/Constrained Plant","Bundle", ...
    "'examples/generated/pcc_flying_plane_contact_n1'");
names = ["q","dq","tip_pose","normal_reaction","contact_feasible","constraint_rcond"];
for index = 1:numel(names)
    addOut(model,names(index),[525 25+42*index 555 45+42*index]);
    add_line(model,"Constrained Plant/"+index,names(index)+"/1");
end
add_line(model,"Arm generalized force/1","Constrained Plant/1");
add_line(model,"Vehicle wrench/1","Constrained Plant/2");
add_line(model,"Tip wrench/1","Constrained Plant/3");
add_line(model,"Constraint acceleration/1","Constrained Plant/4");
save_system(model,outputPath);
close_system(model,0);
end
