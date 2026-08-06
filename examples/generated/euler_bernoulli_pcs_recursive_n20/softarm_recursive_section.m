function y=softarm_recursive_section(section,q,dq,ddqSection,p,v0,w0,a0,alpha0,g0)
%#codegen
switch section
    case 1, y=softarm_recursive_section_template(q(1:2),dq(1:2),ddqSection,p([1,21,41,61,81,101,121,141,161]),v0,w0,a0,alpha0,g0);
    case 2, y=softarm_recursive_section_template(q(3:4),dq(3:4),ddqSection,p([2,22,42,62,82,102,122,142,162]),v0,w0,a0,alpha0,g0);
    case 3, y=softarm_recursive_section_template(q(5:6),dq(5:6),ddqSection,p([3,23,43,63,83,103,123,143,163]),v0,w0,a0,alpha0,g0);
    case 4, y=softarm_recursive_section_template(q(7:8),dq(7:8),ddqSection,p([4,24,44,64,84,104,124,144,164]),v0,w0,a0,alpha0,g0);
    case 5, y=softarm_recursive_section_template(q(9:10),dq(9:10),ddqSection,p([5,25,45,65,85,105,125,145,165]),v0,w0,a0,alpha0,g0);
    case 6, y=softarm_recursive_section_template(q(11:12),dq(11:12),ddqSection,p([6,26,46,66,86,106,126,146,166]),v0,w0,a0,alpha0,g0);
    case 7, y=softarm_recursive_section_template(q(13:14),dq(13:14),ddqSection,p([7,27,47,67,87,107,127,147,167]),v0,w0,a0,alpha0,g0);
    case 8, y=softarm_recursive_section_template(q(15:16),dq(15:16),ddqSection,p([8,28,48,68,88,108,128,148,168]),v0,w0,a0,alpha0,g0);
    case 9, y=softarm_recursive_section_template(q(17:18),dq(17:18),ddqSection,p([9,29,49,69,89,109,129,149,169]),v0,w0,a0,alpha0,g0);
    case 10, y=softarm_recursive_section_template(q(19:20),dq(19:20),ddqSection,p([10,30,50,70,90,110,130,150,170]),v0,w0,a0,alpha0,g0);
    case 11, y=softarm_recursive_section_template(q(21:22),dq(21:22),ddqSection,p([11,31,51,71,91,111,131,151,171]),v0,w0,a0,alpha0,g0);
    case 12, y=softarm_recursive_section_template(q(23:24),dq(23:24),ddqSection,p([12,32,52,72,92,112,132,152,172]),v0,w0,a0,alpha0,g0);
    case 13, y=softarm_recursive_section_template(q(25:26),dq(25:26),ddqSection,p([13,33,53,73,93,113,133,153,173]),v0,w0,a0,alpha0,g0);
    case 14, y=softarm_recursive_section_template(q(27:28),dq(27:28),ddqSection,p([14,34,54,74,94,114,134,154,174]),v0,w0,a0,alpha0,g0);
    case 15, y=softarm_recursive_section_template(q(29:30),dq(29:30),ddqSection,p([15,35,55,75,95,115,135,155,175]),v0,w0,a0,alpha0,g0);
    case 16, y=softarm_recursive_section_template(q(31:32),dq(31:32),ddqSection,p([16,36,56,76,96,116,136,156,176]),v0,w0,a0,alpha0,g0);
    case 17, y=softarm_recursive_section_template(q(33:34),dq(33:34),ddqSection,p([17,37,57,77,97,117,137,157,177]),v0,w0,a0,alpha0,g0);
    case 18, y=softarm_recursive_section_template(q(35:36),dq(35:36),ddqSection,p([18,38,58,78,98,118,138,158,178]),v0,w0,a0,alpha0,g0);
    case 19, y=softarm_recursive_section_template(q(37:38),dq(37:38),ddqSection,p([19,39,59,79,99,119,139,159,179]),v0,w0,a0,alpha0,g0);
    case 20, y=softarm_recursive_section_template(q(39:40),dq(39:40),ddqSection,p([20,40,60,80,100,120,140,160,180]),v0,w0,a0,alpha0,g0);
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
end
