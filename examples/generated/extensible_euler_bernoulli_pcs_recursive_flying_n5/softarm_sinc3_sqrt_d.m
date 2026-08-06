function y = softarm_sinc3_sqrt_d(z)
%#codegen
if abs(z)<1e-8, y=-1/120+z/2520-z^2/120960+z^3/9979200; else, a=softarm_sinc_sqrt(z); y=(a-1-z*softarm_sinc_sqrt_d(z))/z^2; end
end
