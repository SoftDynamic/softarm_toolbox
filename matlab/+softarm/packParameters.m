function p = packParameters(manifest, values)
%PACKPARAMETERS Pack named overrides using manifest order.
arguments
    manifest (1,1) struct
    values (1,1) struct = struct()
end
p = reshape([manifest.parameters.default], [], 1);
for index = 1:numel(manifest.parameters)
    name = manifest.parameters(index).name;
    if isfield(values, name)
        value = values.(name);
        validateattributes(value, {'numeric'}, {'real', 'finite', 'scalar'}, ...
            "softarm.packParameters", name);
        p(index) = value;
    end
end
unknown = setdiff(fieldnames(values), {manifest.parameters.name});
assert(isempty(unknown), "softarm:UnknownParameter", ...
    "Unknown model parameter(s): %s", strjoin(unknown, ", "));
end
