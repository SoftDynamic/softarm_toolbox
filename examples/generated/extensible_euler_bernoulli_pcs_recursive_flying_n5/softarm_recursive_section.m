function y=softarm_recursive_section(section,q,dq,ddqSection,p,v0,w0,a0,alpha0,g0)
%#codegen
switch section
    case 1, y=softarm_recursive_section_template(q(7:9),dq(7:9),ddqSection,p([1,6,11,16,21,26,31,36,41,46,51]),v0,w0,a0,alpha0,g0);
    case 2, y=softarm_recursive_section_template(q(10:12),dq(10:12),ddqSection,p([2,7,12,17,22,27,32,37,42,47,52]),v0,w0,a0,alpha0,g0);
    case 3, y=softarm_recursive_section_template(q(13:15),dq(13:15),ddqSection,p([3,8,13,18,23,28,33,38,43,48,53]),v0,w0,a0,alpha0,g0);
    case 4, y=softarm_recursive_section_template(q(16:18),dq(16:18),ddqSection,p([4,9,14,19,24,29,34,39,44,49,54]),v0,w0,a0,alpha0,g0);
    case 5, y=softarm_recursive_section_template(q(19:21),dq(19:21),ddqSection,p([5,10,15,20,25,30,35,40,45,50,55]),v0,w0,a0,alpha0,g0);
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
end
