function y=softarm_recursive_end(section,q,dq,p)
%#codegen
switch section
    case 1, y=softarm_recursive_end_template(q(1:2),dq(1:2),p([1,21,41,61,81,101,121,141,161]));
    case 2, y=softarm_recursive_end_template(q(3:4),dq(3:4),p([2,22,42,62,82,102,122,142,162]));
    case 3, y=softarm_recursive_end_template(q(5:6),dq(5:6),p([3,23,43,63,83,103,123,143,163]));
    case 4, y=softarm_recursive_end_template(q(7:8),dq(7:8),p([4,24,44,64,84,104,124,144,164]));
    case 5, y=softarm_recursive_end_template(q(9:10),dq(9:10),p([5,25,45,65,85,105,125,145,165]));
    case 6, y=softarm_recursive_end_template(q(11:12),dq(11:12),p([6,26,46,66,86,106,126,146,166]));
    case 7, y=softarm_recursive_end_template(q(13:14),dq(13:14),p([7,27,47,67,87,107,127,147,167]));
    case 8, y=softarm_recursive_end_template(q(15:16),dq(15:16),p([8,28,48,68,88,108,128,148,168]));
    case 9, y=softarm_recursive_end_template(q(17:18),dq(17:18),p([9,29,49,69,89,109,129,149,169]));
    case 10, y=softarm_recursive_end_template(q(19:20),dq(19:20),p([10,30,50,70,90,110,130,150,170]));
    case 11, y=softarm_recursive_end_template(q(21:22),dq(21:22),p([11,31,51,71,91,111,131,151,171]));
    case 12, y=softarm_recursive_end_template(q(23:24),dq(23:24),p([12,32,52,72,92,112,132,152,172]));
    case 13, y=softarm_recursive_end_template(q(25:26),dq(25:26),p([13,33,53,73,93,113,133,153,173]));
    case 14, y=softarm_recursive_end_template(q(27:28),dq(27:28),p([14,34,54,74,94,114,134,154,174]));
    case 15, y=softarm_recursive_end_template(q(29:30),dq(29:30),p([15,35,55,75,95,115,135,155,175]));
    case 16, y=softarm_recursive_end_template(q(31:32),dq(31:32),p([16,36,56,76,96,116,136,156,176]));
    case 17, y=softarm_recursive_end_template(q(33:34),dq(33:34),p([17,37,57,77,97,117,137,157,177]));
    case 18, y=softarm_recursive_end_template(q(35:36),dq(35:36),p([18,38,58,78,98,118,138,158,178]));
    case 19, y=softarm_recursive_end_template(q(37:38),dq(37:38),p([19,39,59,79,99,119,139,159,179]));
    case 20, y=softarm_recursive_end_template(q(39:40),dq(39:40),p([20,40,60,80,100,120,140,160,180]));
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
end
