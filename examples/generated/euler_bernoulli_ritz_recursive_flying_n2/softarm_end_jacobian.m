function J=softarm_end_jacobian(q,p)
%#codegen
[B,dB,d2B]=softarm_affine_base_jet(q,p); qa=q(7:end);
A0=softarm_affine_mount(p); dA=zeros(4,4,4);

for section=1:2
    first=(section-1)*2+1; last=first+2-1; [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last);
end
[~,Jv,~,Jw,~,~]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa); J=[Jv;Jw];
end
