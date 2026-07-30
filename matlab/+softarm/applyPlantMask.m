function plant = applyPlantMask(owner, requirement)
%APPLYPLANTMASK Resolve and apply a softarm_plant System Mask argument.
arguments
    owner (1,1) string
    requirement (1,1) string = "plant"
end

if string(get_param(owner,"Type")) == "block_diagram"
    modelWorkspace = get_param(owner,"ModelWorkspace");
    bundle = string(getVariable(modelWorkspace,"Bundle"));
else
    expression = char(get_param(owner,"Bundle"));
    try
        bundle = string(evalin("base",expression));
    catch cause
        if contains(expression,"/") || contains(expression,"\")
            bundle = string(expression);
        else
            rethrow(cause);
        end
    end
end
plant = softarm.configureBundle(owner,bundle,requirement);
end
