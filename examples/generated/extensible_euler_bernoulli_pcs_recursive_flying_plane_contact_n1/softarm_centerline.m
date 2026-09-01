function P=softarm_centerline(q,p,xi)
%#codegen
xi=double(xi(:).'); assert(~isempty(xi)&&all(isfinite(xi))&&all(xi>=0)&&all(xi<=1));
P=zeros(3,1*numel(xi));
for sample=1:numel(xi)
    points=softarm_centerline_at(q,p,xi(sample));
    for section=1:1, P(:,(section-1)*numel(xi)+sample)=points(:,section); end
end
end
