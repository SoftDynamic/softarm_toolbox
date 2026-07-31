function y = softarm_sinc_sqrt_d(z)
%#codegen
if abs(z)<1e-8, y=-1/6+z/60-z^2/1680+z^3/90720; else, s=sqrt(z); y=(s*cos(s)-sin(s))/(2*s^3); end
end
