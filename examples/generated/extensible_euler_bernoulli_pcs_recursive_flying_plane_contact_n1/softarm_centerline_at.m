function P=softarm_centerline_at(q,p,xi)
%#codegen
P=zeros(3,1); current=softarm_mount_transform(q,p);
for section=1:1
    material=current*softarm_recursive_pose(section,q,p,xi); P(:,section)=material(1:3,4);
    current=current*softarm_recursive_pose(section,q,p,1);
end
end
