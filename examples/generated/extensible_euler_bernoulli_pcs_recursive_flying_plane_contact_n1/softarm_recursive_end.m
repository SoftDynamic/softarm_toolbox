function y=softarm_recursive_end(section,q,dq,p)
%#codegen
switch section
    case 1, y=softarm_recursive_end_template(q(7:9),dq(7:9),p([1,2,3,4,5,6,7,8,9,10,11]));
    otherwise, error('softarm:InvalidSection','Invalid section index.');
end
end
