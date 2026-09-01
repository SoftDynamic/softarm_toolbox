function [B,dB,d2B]=softarm_affine_base_jet(q,p)
%#codegen
y=softarm_affine_base_jet_raw(q,p); cursor=1;
B=reshape(y(cursor:cursor+15),4,4); cursor=cursor+16;
dB=zeros(4,4,6); d2B=zeros(4,4,6,6);
for first=1:6
    dB(:,:,first)=reshape(y(cursor:cursor+15),4,4); cursor=cursor+16;
end
for first=1:6
    for second=1:6
        d2B(:,:,first,second)=reshape(y(cursor:cursor+15),4,4); cursor=cursor+16;
    end
end
end
