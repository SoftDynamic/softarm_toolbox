function info = softarm_actuator_info()
%SOFTARM_ACTUATOR_INFO Ordered runtime actuator information.
info=struct('names',["bend_x","bend_y"],'kinds',["signed","signed"],'count',2,'acceleration',"strict");
end
