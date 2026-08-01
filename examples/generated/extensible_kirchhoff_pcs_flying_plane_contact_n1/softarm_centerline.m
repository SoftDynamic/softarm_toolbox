function P = softarm_centerline(q,p,xi)
%SOFTARM_CENTERLINE Sample model-derived section centerlines.
xi=double(xi(:).');
assert(~isempty(xi)&&all(isfinite(xi))&&all(xi>=0)&&all(xi<=1),'softarm:InvalidMaterialCoordinate','xi must be finite and lie in [0,1].');
P=zeros(3,1*numel(xi));
for sample=1:numel(xi)
    sectionPoints=softarm_centerline_at(q,p,xi(sample));
    for section=1:1
        P(:,(section-1)*numel(xi)+sample)=sectionPoints(:,section);
    end
end
end
