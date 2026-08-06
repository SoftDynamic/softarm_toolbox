function y = softarm_sinc3_sqrt(z)
%#codegen
if abs(z)<1e-8, y=1/6-z/120+z^2/5040-z^3/362880+z^4/39916800; else, y=(1-softarm_sinc_sqrt(z))/z; end
end
