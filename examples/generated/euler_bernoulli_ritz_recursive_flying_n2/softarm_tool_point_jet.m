function [position,J,gamma,velocity]=softarm_tool_point_jet(q,dq,toolOffset,p)
%#codegen
[B,dB,d2B]=softarm_affine_base_jet(q,p); qa=q(7:end);
A0=softarm_affine_mount(p); dA=zeros(4,4,4);

for section=1:2
    first=(section-1)*2+1; last=first+2-1; [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last);
end
tool=[eye(3),toolOffset(:);0 0 0 1]; A0=A0*tool;
for arm=1:4, dA(:,:,arm)=dA(:,:,arm)*tool; end
[T,J,dJ,~,~,~]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa); Jdot=zeros(3,10); for direction=1:10, Jdot=Jdot+dJ(:,:,direction)*dq(direction); end
position=T(1:3,4); velocity=J*dq; gamma=Jdot*dq;
end
