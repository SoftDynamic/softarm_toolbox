function info = softarm_actuator_info()
%SOFTARM_ACTUATOR_INFO Ordered runtime actuator information.
info=struct('names',["t1","t2","t3"],'kinds',["unilateral","unilateral","unilateral"],'count',3,'acceleration',"strict");
end
