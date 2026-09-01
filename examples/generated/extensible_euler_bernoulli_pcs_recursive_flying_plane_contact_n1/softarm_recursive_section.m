function y=softarm_recursive_section(section,q,dq,ddqSection,p,v0,w0,a0,alpha0,g0)
%#codegen
switch section
    case 1, y=softarm_recursive_section_template(q(7:9),dq(7:9),ddqSection,p([1,2,3,4,5,6,7,8,9,10,11]),v0,w0,a0,alpha0,g0);
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
end
