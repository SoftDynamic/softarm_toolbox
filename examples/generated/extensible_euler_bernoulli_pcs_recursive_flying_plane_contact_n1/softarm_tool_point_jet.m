function [position,J,gamma,velocity]=softarm_tool_point_jet(q,dq,toolOffset,p)
%SOFTARM_TOOL_POINT_JET World-frame point kinematics and acceleration bias.
%#codegen
q=q(:); dq=dq(:); toolOffset=toolOffset(:); assert(numel(toolOffset)==3);
base=softarm_recursive_base(q,dq,zeros(6,1),p);
v=base(1:3); w=base(4:6); a=base(7:9); alpha=base(10:12);
T=softarm_mount_transform(q,p);
for section=1:1
    first=6+(section-1)*3+1; last=first+3-1;
    jet=softarm_recursive_end(section,q,dq,p); cursor=1;
    rotation=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    offset=jet(cursor:cursor+2); cursor=cursor+3;
    Jv=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    Jw=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    Jvd=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    Jwd=reshape(jet(cursor:cursor+8),3,3);
    relativeV=Jv*dq(first:last); relativeW=Jw*dq(first:last);
    vParent=v+cross(w,offset)+relativeV; wParent=w+relativeW;
    aParent=a+cross(alpha,offset)+cross(w,cross(w,offset))+2*cross(w,relativeV)+Jvd*dq(first:last);
    alphaParent=alpha+cross(w,relativeW)+Jwd*dq(first:last);
    v=rotation.'*vParent; w=rotation.'*wParent;
    a=rotation.'*aParent; alpha=rotation.'*alphaParent;
    T=T*[rotation offset;0 0 0 1];
end
R=T(1:3,1:3); worldOffset=R*toolOffset;
position=T(1:3,4)+worldOffset;
endJ=softarm_end_jacobian(q,p);
J=endJ(1:3,:)-softarm_skew(worldOffset)*endJ(4:6,:);
velocity=J*dq;
worldW=R*w; worldA=R*a; worldAlpha=R*alpha;
gamma=worldA+cross(worldAlpha,worldOffset)+cross(worldW,cross(worldW,worldOffset));
end
