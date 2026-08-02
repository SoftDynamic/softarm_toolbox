function [c,s] = softarm_affine_moment_parts(n,c0,c1,xi)
%#codegen
assert(n>=0&&n==floor(n));
if xi==0, c=0; s=0; return; end
previousPreviousReal=0; previousPreviousImag=0;
previousReal=1; previousImag=0;
c=xi^(n+1)/(n+1); s=0; smallTerms=0;
for degree=1:256
    realCoefficient=(-c0*previousImag-c1*previousPreviousImag)/degree;
    imagCoefficient=(c0*previousReal+c1*previousPreviousReal)/degree;
    scale=xi^(degree+n+1)/(degree+n+1);
    realTerm=realCoefficient*scale; imagTerm=imagCoefficient*scale;
    c=c+realTerm; s=s+imagTerm;
    if hypot(realTerm,imagTerm)<=2e-16*max(1,hypot(c,s))
        smallTerms=smallTerms+1;
        if smallTerms>=4, break; end
    else
        smallTerms=0;
    end
    previousPreviousReal=previousReal; previousPreviousImag=previousImag;
    previousReal=realCoefficient; previousImag=imagCoefficient;
end
end
