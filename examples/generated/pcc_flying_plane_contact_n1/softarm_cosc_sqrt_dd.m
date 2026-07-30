function y = softarm_cosc_sqrt_dd(z)
%#codegen
if abs(z)<1e-8, y=1/360-z/6720+z^2/302400; else, s=sqrt(z); y=(z*cos(s)-5*s*sin(s)+8-8*cos(s))/(4*z^3); end
end
