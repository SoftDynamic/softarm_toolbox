function phi=softarm_constraint_value(q,p)
[phi,~,~,~]=softarm_constraint_terms(q,zeros(9,1),p);
end
