function A=softarm_constraint_jacobian(q,p)
[~,A,~,~]=softarm_constraint_terms(q,zeros(9,1),p);
end
