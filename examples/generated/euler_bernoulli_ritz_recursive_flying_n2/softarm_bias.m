function h=softarm_bias(q,dq,p)
%#codegen
h=softarm_inverse_dynamics(q,dq,zeros(10,1),p);
end
