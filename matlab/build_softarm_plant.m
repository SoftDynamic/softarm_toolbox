function build_softarm_plant(outputPath)
%BUILD_SOFTARM_PLANT Rebuild the repository's generic Simulink plant.
arguments
    outputPath (1,1) string = fullfile(fileparts(fileparts(mfilename("fullpath"))), "softarm_plant.slx")
end
model = "softarm_plant";
if bdIsLoaded(model), close_system(model, 0); end
new_system(model);
defaultBundle = "examples/generated/pcc_lumped_n2";
modelWorkspace = get_param(model,"ModelWorkspace");
assignin(modelWorkspace,"Bundle",char(defaultBundle));
set_param(model,"ParameterArgumentNames","Bundle");
set_param(model, "Solver", "ode15s", "StopTime", "5");
set_param(model,"PreLoadFcn","softarm_root=fileparts(get_param(bdroot,'FileName'));addpath(softarm_root);addpath(fullfile(softarm_root,'matlab'));");
mask = Simulink.Mask.create(model);
mask.Description = "Select a generated SoftArm bundle. Changing the bundle recompiles dimensions and parameters.";
bundleParameter = mask.getParameter("Bundle");
bundleParameter.Prompt = "Bundle directory";
bundleParameter.Callback = "if string(get_param(gcb,'Type'))~=""block_diagram"",softarm.applyPlantMask(string(gcb),'plant');end;";

add_block("simulink/Sources/In1", model+"/tau_arm", "Position", [30 55 60 75], "PortDimensions", "softarm_narm");
add_block("simulink/Sources/In1", model+"/vehicle_wrench", "Position", [30 105 60 125], "PortDimensions", "6");
add_block("simulink/Sources/In1", model+"/tip_wrench", "Position", [30 155 60 175], "PortDimensions", "6");
add_block("simulink/Sources/Constant", model+"/parameters", "Position", [30 215 90 245], "Value", "softarm_p");
add_block("simulink/Continuous/Integrator", model+"/state", "Position", [360 75 390 105], "InitialCondition", "softarm_x0");

add_block("simulink/Signal Routing/Demux", model+"/State split", "Position", [440 55 445 145], ...
    "Outputs", "[softarm_nq softarm_nq]");

add_block("simulink/User-Defined Functions/MATLAB Function", model+"/State RHS", "Position", [170 65 300 175]);
chart = find(sfroot, "-isa", "Stateflow.EMChart", "Path", model+"/State RHS");
chart.Script = sprintf("function dx = fcn(x,tauArm,wVehicle,wTip,p)\n%%#codegen\ndx = softarm_state_rhs(x,tauArm,wVehicle,wTip,p);\nend\n");

add_block("simulink/User-Defined Functions/MATLAB Function", model+"/Tip pose", "Position", [540 170 650 230]);
chart = find(sfroot, "-isa", "Stateflow.EMChart", "Path", model+"/Tip pose");
chart.Script = sprintf("function pose = fcn(q,p)\n%%#codegen\nH=softarm_kinematics(q,p); pose=reshape(H(:,:,end),16,1);\nend\n");

add_block("simulink/User-Defined Functions/MATLAB Function", model+"/Diagnostics", "Position", [540 255 650 315]);
chart = find(sfroot, "-isa", "Stateflow.EMChart", "Path", model+"/Diagnostics");
chart.Script = sprintf("function diagnostic = fcn(q,p)\n%%#codegen\ndiagnostic=rcond(softarm_mass(q,p));\nend\n");

add_block("simulink/Sinks/Out1", model+"/q_out", "Position", [720 55 750 75]);
add_block("simulink/Sinks/Out1", model+"/dq_out", "Position", [720 110 750 130]);
add_block("simulink/Sinks/Out1", model+"/tip_pose", "Position", [720 190 750 210]);
add_block("simulink/Sinks/Out1", model+"/diagnostic", "Position", [720 275 750 295]);

add_line(model, "state/1", "State split/1");
add_line(model, "state/1", "State RHS/1"); add_line(model, "tau_arm/1", "State RHS/2");
add_line(model, "vehicle_wrench/1", "State RHS/3"); add_line(model, "tip_wrench/1", "State RHS/4");
add_line(model, "parameters/1", "State RHS/5");
add_line(model, "State RHS/1", "state/1");
add_line(model, "State split/1", "q_out/1"); add_line(model, "State split/2", "dq_out/1");
add_line(model, "State split/1", "Tip pose/1"); add_line(model, "parameters/1", "Tip pose/2");
add_line(model, "Tip pose/1", "tip_pose/1");
add_line(model, "State split/1", "Diagnostics/1"); add_line(model, "parameters/1", "Diagnostics/2");
add_line(model, "Diagnostics/1", "diagnostic/1");

save_system(model, outputPath);
close_system(model, 0);
end
