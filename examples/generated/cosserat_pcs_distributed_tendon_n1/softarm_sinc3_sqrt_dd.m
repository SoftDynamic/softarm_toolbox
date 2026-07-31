function y = softarm_sinc3_sqrt_dd(z)
%#codegen
if abs(z)<1e-8, y=1/2520-z/60480+z^2/3326400; else, a=softarm_sinc_sqrt(z); y=(2-2*a+2*z*softarm_sinc_sqrt_d(z)-z^2*softarm_sinc_sqrt_dd(z))/z^3; end
end
