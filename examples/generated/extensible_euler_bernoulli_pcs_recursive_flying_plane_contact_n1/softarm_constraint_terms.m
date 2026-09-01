function [phi,A,gamma,G]=softarm_constraint_terms(q,dq,p)
%#codegen
offset=p([24 25 26]);
[pointPosition,pointJacobian,pointVelocityBias,pointVelocity]=softarm_tool_point_jet(q,dq,offset,p);
y=softarm_point_constraint_kernel(pointPosition,pointJacobian,pointVelocityBias,pointVelocity,p);
cursor=1; phi=y(cursor); cursor=cursor+1;
A=reshape(y(cursor:cursor+8),1,9); cursor=cursor+9;
gamma=y(cursor); cursor=cursor+1;
G=reshape(y(cursor:20),9,1);
end
