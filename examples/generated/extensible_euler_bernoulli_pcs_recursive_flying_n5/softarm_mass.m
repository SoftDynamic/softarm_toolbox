function M=softarm_mass(q,p)
%#codegen
n=21; z=zeros(n,1); reference=softarm_inverse_dynamics(q,z,z,p); M=zeros(n,n);
for column=1:n
    acceleration=zeros(n,1); acceleration(column)=1;
    M(:,column)=softarm_inverse_dynamics(q,z,acceleration,p)-reference;
end
M=(M+M.')/2;
end
