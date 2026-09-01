function h=softarm_bias(q,dq,p)
%#codegen
h=softarm_inverse_dynamics(q,dq,zeros(9,1),p);
end
