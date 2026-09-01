function [F0,dF]=softarm_legacy_affine(section,p,xi)
%#codegen
switch section
    case 1, y=softarm_legacy_affine_template(p([1,3,5,7,9,11,13,15,17]),xi);
    case 2, y=softarm_legacy_affine_template(p([2,4,6,8,10,12,14,16,18]),xi);
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
F0=reshape(y(1:16),4,4); dF=reshape(y(17:end),4,4,2);
end
