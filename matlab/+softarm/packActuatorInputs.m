function input = packActuatorInputs(actuation, values)
%PACKACTUATORINPUTS Pack named actuator channel values in generated order.
arguments
    actuation (1,1) struct
    values (1,1) struct = struct()
end
names = cellstr(actuation.names);
input = zeros(actuation.count,1);
for index = 1:actuation.count
    name = names{index};
    if isfield(values,name)
        value = values.(name);
        validateattributes(value,{'numeric'},{'real','finite','scalar'}, ...
            "softarm.packActuatorInputs",name);
        input(index) = value;
    end
end
unknown = setdiff(fieldnames(values),names);
assert(isempty(unknown),"softarm:UnknownActuatorChannel", ...
    "Unknown actuator channel(s): %s",strjoin(unknown,", "));
end
