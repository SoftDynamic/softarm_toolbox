function y = softarm_cosc_sqrt(z)
%#codegen
if abs(z)<1e-8, y=1/2-z/24+z^2/720-z^3/40320+z^4/3628800; else, s=sqrt(z); y=(1-cos(s))/z; end
end
