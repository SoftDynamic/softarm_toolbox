function H=softarm_kinematics(q,p)
%#codegen
H=zeros(4,4,1); current=softarm_mount_transform(q,p);
for section=1:1
    current=current*softarm_recursive_pose(section,q,p,1); H(:,:,section)=current;
end
end
