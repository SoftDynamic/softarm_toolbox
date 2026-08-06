function J=softarm_end_jacobian(q,p)
%#codegen
J=softarm_mount_jacobian(q,p); T=softarm_mount_transform(q,p); zeroDq=zeros(21,1);
for section=1:5
    jet=softarm_recursive_end(section,q,zeroDq,p); cursor=1;
    rotation=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    position=jet(cursor:cursor+2); cursor=cursor+3;
    localJv=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    localJw=reshape(jet(cursor:cursor+8),3,3);
    worldOffset=T(1:3,1:3)*position;
    J(1:3,:)=J(1:3,:)-softarm_skew(worldOffset)*J(4:6,:);
    first=6+(section-1)*3+1; last=first+3-1;
    J(1:3,first:last)=J(1:3,first:last)+T(1:3,1:3)*localJv;
    J(4:6,first:last)=J(4:6,first:last)+T(1:3,1:3)*localJw;
    T=T*[rotation position;0 0 0 1];
end
end
