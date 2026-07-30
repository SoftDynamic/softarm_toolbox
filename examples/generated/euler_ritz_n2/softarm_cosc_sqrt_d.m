function y = softarm_cosc_sqrt_d(z)
%#codegen
if abs(z)<1e-8, y=-1/24+z/360-z^2/13440+z^3/907200; else, s=sqrt(z); y=(s*sin(s)-2*(1-cos(s)))/(2*z^2); end
end
