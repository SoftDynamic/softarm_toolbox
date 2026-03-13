function S = skew3(v)
arguments (Input)
    v (3,1)
end
arguments (Output)
    S (3,3)
end

S = [0      -v(3)       v(2);
     v(3)   0           -v(1);
     -v(2)  v(1)        0];
end