function y = softarm_sinc_sqrt_dd(z)
%#codegen
if abs(z)<1e-8, y=1/60-z/840+z^2/30240; else, s=sqrt(z); y=((3-z)*sin(s)-3*s*cos(s))/(4*s^5); end
end
