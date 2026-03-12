function v = vex3(S)
arguments (Input)
    S (3,3)
end
arguments (Output)
    v (3,1)
end

% assert(S == -S', 'Input is not a skew-symmetric matrix');
v = [S(3,2); S(1,3); S(2,1)];
end