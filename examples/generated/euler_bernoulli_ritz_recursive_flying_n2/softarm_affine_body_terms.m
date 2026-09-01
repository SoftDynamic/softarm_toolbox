function [M,Mdot,gradK,gravityTau]=softarm_affine_body_terms(B,dB,d2B,A0,dA,qa,dq,mass,inertia,gravity)
%#codegen
[~,Jv,dJv,Jw,dJw,R0]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa);
n=numel(dq); worldInertia=R0*inertia*R0.'; M=mass*(Jv.'*Jv)+Jw.'*worldInertia*Jw;
Mdot=zeros(n); gradK=zeros(n,1); nb=size(dB,3);
for direction=1:n
    dWorldInertia=zeros(3);
    if direction<=nb
        dR=dB(1:3,1:3,direction)*A0(1:3,1:3);
        dWorldInertia=dR*inertia*R0.'+R0*inertia*dR.';
    end
    dM=mass*(dJv(:,:,direction).'*Jv+Jv.'*dJv(:,:,direction))+dJw(:,:,direction).'*worldInertia*Jw+Jw.'*dWorldInertia*Jw+Jw.'*worldInertia*dJw(:,:,direction);
    Mdot=Mdot+dM*dq(direction); gradK(direction)=0.5*dq.'*dM*dq;
end
gravityTau=-mass*gravity*Jv(3,:).';
end
