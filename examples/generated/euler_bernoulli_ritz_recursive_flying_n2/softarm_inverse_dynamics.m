function tau=softarm_inverse_dynamics(q,dq,ddq,p)
%#codegen
[M,h]=softarm_affine_components(q,dq,p); tau=M*ddq(:)+h;
end
