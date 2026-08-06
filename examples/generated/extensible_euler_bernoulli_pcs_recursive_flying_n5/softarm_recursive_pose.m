function H=softarm_recursive_pose(section,q,p,xi)
%#codegen
switch section
    case 1, H=softarm_recursive_pose_template(q(7:9),p([1,6,11,16,21,26,31,36,41,46,51]),xi);
    case 2, H=softarm_recursive_pose_template(q(10:12),p([2,7,12,17,22,27,32,37,42,47,52]),xi);
    case 3, H=softarm_recursive_pose_template(q(13:15),p([3,8,13,18,23,28,33,38,43,48,53]),xi);
    case 4, H=softarm_recursive_pose_template(q(16:18),p([4,9,14,19,24,29,34,39,44,49,54]),xi);
    case 5, H=softarm_recursive_pose_template(q(19:21),p([5,10,15,20,25,30,35,40,45,50,55]),xi);
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
end
