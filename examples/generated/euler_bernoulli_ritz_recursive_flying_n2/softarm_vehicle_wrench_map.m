function Bv=softarm_vehicle_wrench_map(q,p)
%#codegen
[B,dB,d2B]=softarm_affine_base_jet(q,p); [~,Jv,~,Jw,~,~]=softarm_affine_kinematic_jet(B,dB,d2B,eye(4),zeros(4,4,4),q(7:end)); J=[Jv;Jw]; R=B(1:3,1:3); Bv=J.'*blkdiag(R,R);
end
