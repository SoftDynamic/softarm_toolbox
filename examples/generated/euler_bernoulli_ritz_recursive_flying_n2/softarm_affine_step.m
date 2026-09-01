function [next0,dNext]=softarm_affine_step(parent0,dParent,F0,dF,first,last)
%#codegen
n=size(dParent,3); dNext=zeros(4,4,n); next0=parent0*F0;
for coordinate=1:n, dNext(:,:,coordinate)=dParent(:,:,coordinate)*F0; end
for local=1:size(dF,3), dNext(:,:,first+local-1)=parent0*dF(:,:,local); end
end
