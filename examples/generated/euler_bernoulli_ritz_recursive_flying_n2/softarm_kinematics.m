function H=softarm_kinematics(q,p)
%#codegen
[B,dB,d2B]=softarm_affine_base_jet(q,p); qa=q(7:end);
A0=softarm_affine_mount(p); dA=zeros(4,4,4);
H=zeros(4,4,2);
for section=1:2
    first=(section-1)*2+1; last=first+2-1; [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last); A=A0;
    for arm=1:4, A=A+dA(:,:,arm)*qa(arm); end
    H(:,:,section)=B*A;
end
end
