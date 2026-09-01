function H=softarm_recursive_pose(section,q,p,xi)
%#codegen
switch section
    case 1, H=softarm_recursive_pose_template(q(7:9),p([1,2,3,4,5,6,7,8,9,10,11]),xi);
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
end
