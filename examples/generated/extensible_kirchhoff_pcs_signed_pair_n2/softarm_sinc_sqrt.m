function y = softarm_sinc_sqrt(z)
%#codegen
if abs(z)<1e-8, y=1-z/6+z^2/120-z^3/5040+z^4/362880; else, s=sqrt(z); y=sin(s)/s; end
end
