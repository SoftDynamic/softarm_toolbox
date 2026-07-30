function info = softarm_actuator_info()
%SOFTARM_ACTUATOR_INFO Ordered runtime actuator information.
info=struct('names',["pair_x","pair_y"],'kinds',["signed","signed"],'count',2,'acceleration',"strict");
end
