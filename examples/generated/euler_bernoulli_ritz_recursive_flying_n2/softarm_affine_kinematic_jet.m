function [T,Jv,dJv,Jw,dJw,R0]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa)
%#codegen
nb=size(dB,3); na=size(dA,3); n=nb+na; A=A0;
for arm=1:na, A=A+dA(:,:,arm)*qa(arm); end
T=B*A; R0=B(1:3,1:3)*A0(1:3,1:3); Jv=zeros(3,n); dJv=zeros(3,n,n);
D=zeros(3,3,n); dD=zeros(3,3,n,n); dR0=zeros(3,3,n);
for column=1:nb
    derivative=dB(:,:,column)*A; Jv(:,column)=derivative(1:3,4);
    D(:,:,column)=dB(1:3,1:3,column)*A0(1:3,1:3);
    dR0(:,:,column)=D(:,:,column);
end
for arm=1:na
    column=nb+arm; derivative=B*dA(:,:,arm); Jv(:,column)=derivative(1:3,4);
    D(:,:,column)=B(1:3,1:3)*dA(1:3,1:3,arm);
end
for direction=1:nb
    for column=1:nb
        derivative=d2B(:,:,column,direction)*A;
        dJv(:,column,direction)=derivative(1:3,4);
        dD(:,:,column,direction)=d2B(1:3,1:3,column,direction)*A0(1:3,1:3);
    end
    for arm=1:na
        column=nb+arm; derivative=dB(:,:,direction)*dA(:,:,arm);
        dJv(:,column,direction)=derivative(1:3,4);
        dD(:,:,column,direction)=dB(1:3,1:3,direction)*dA(1:3,1:3,arm);
    end
end
for direction=1:na
    for column=1:nb
        derivative=dB(:,:,column)*dA(:,:,direction);
        dJv(:,column,nb+direction)=derivative(1:3,4);
    end
end
Jw=zeros(3,n); dJw=zeros(3,n,n);
for column=1:n
    rate=D(:,:,column)*R0.'; skew=(rate-rate.')/2;
    Jw(:,column)=[skew(3,2);skew(1,3);skew(2,1)];
    for direction=1:nb
        rateDerivative=dD(:,:,column,direction)*R0.'+D(:,:,column)*dR0(:,:,direction).';
        skew=(rateDerivative-rateDerivative.')/2;
        dJw(:,column,direction)=[skew(3,2);skew(1,3);skew(2,1)];
    end
end
end
