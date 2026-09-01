function [M,h]=softarm_affine_components(q,dq,p)
%#codegen
q=q(:); dq=dq(:); qa=q(7:end);
[B,dB,d2B]=softarm_affine_base_jet(q,p); A0=softarm_affine_mount(p); dA=zeros(4,4,4);
nodes=[0.0694318442029737123880267555536;0.330009478207571867598667120448;0.669990521792428132401332879552;0.930568155797026287611973244446]; weights=[0.17392742256872693;0.32607257743127307;0.32607257743127307;0.17392742256872693];
massIndex=[3;4]; ixxIndex=[5;6]; iyyIndex=[7;8]; izzIndex=[9;10];
M=zeros(10); Mdot=zeros(10); gradK=zeros(10,1); gravityTau=zeros(10,1);

vehicleInertia=diag([p(25),p(26),p(27)]);
[bodyM,bodyMdot,bodyGrad,bodyGravity]=softarm_affine_body_terms(B,dB,d2B,eye(4),zeros(4,4,4),qa,dq,p(24),vehicleInertia,p(19));
M=M+bodyM; Mdot=Mdot+bodyMdot; gradK=gradK+bodyGrad; gravityTau=gravityTau+bodyGravity;

for section=1:2
    first=(section-1)*2+1; last=first+2-1;
    inertia=diag([p(ixxIndex(section)),p(iyyIndex(section)),p(izzIndex(section))]);
    for sample=1:numel(nodes)
        [F0,dF]=softarm_legacy_affine(section,p,nodes(sample));
        [material0,dMaterial]=softarm_affine_step(A0,dA,F0,dF,first,last);
        [bodyM,bodyMdot,bodyGrad,bodyGravity]=softarm_affine_body_terms(B,dB,d2B,material0,dMaterial,qa,dq,p(massIndex(section)),inertia,p(19));
        M=M+weights(sample)*bodyM; Mdot=Mdot+weights(sample)*bodyMdot;
        gradK=gradK+weights(sample)*bodyGrad; gravityTau=gravityTau+weights(sample)*bodyGravity;
    end
    [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last);
end
tipInertia=diag([p(21),p(22),p(23)]);
[bodyM,bodyMdot,bodyGrad,bodyGravity]=softarm_affine_body_terms(B,dB,d2B,A0,dA,qa,dq,p(20),tipInertia,p(19));
M=M+bodyM; Mdot=Mdot+bodyMdot; gradK=gradK+bodyGrad; gravityTau=gravityTau+bodyGravity;
M=(M+M.')/2; h=Mdot*dq-gradK+gravityTau+softarm_recursive_internal_force(q,dq,p);
end
