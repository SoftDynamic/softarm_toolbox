function tau=softarm_inverse_dynamics(q,dq,ddq,p)
%SOFTARM_INVERSE_DYNAMICS Recursive section inverse dynamics.
%#codegen
assert(numel(q)==21&&numel(dq)==21&&numel(ddq)==21);
q=q(:); dq=dq(:); ddq=ddq(:);
base=softarm_recursive_base(q,dq,ddq(1:6),p);
assert(numel(base)==57);
v=base(1:3); w=base(4:6); a=base(7:9); alpha=base(10:12); g=base(13:15);
baseTau=base(16:15+6);
mountMap=reshape(base(16+6:end),6,6);
Rend=zeros(3,3,5); rend=zeros(3,5);
Jv=zeros(3,3,5); Jw=zeros(3,3,5);
ownWrench=zeros(6,5); ownTau=zeros(3,5);
for section=1:5
    first=6+(section-1)*3+1; last=first+3-1;
    jet=softarm_recursive_end(section,q,dq,p);
    cursor=1; Rend(:,:,section)=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    rend(:,section)=jet(cursor:cursor+2); cursor=cursor+3;
    Jv(:,:,section)=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    Jw(:,:,section)=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    Jvd=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    Jwd=reshape(jet(cursor:cursor+8),3,3);
    local=softarm_recursive_section(section,q,dq,ddq(first:last),p,v,w,a,alpha,g);
    ownWrench(:,section)=local(1:6); ownTau(:,section)=local(7:end);
    relativeV=Jv(:,:,section)*dq(first:last); relativeW=Jw(:,:,section)*dq(first:last);
    vParent=v+cross(w,rend(:,section))+relativeV;
    wParent=w+relativeW;
    aParent=a+cross(alpha,rend(:,section))+cross(w,cross(w,rend(:,section)))+2*cross(w,relativeV)+Jv(:,:,section)*ddq(first:last)+Jvd*dq(first:last);
    alphaParent=alpha+cross(w,relativeW)+Jw(:,:,section)*ddq(first:last)+Jwd*dq(first:last);
    rotation=Rend(:,:,section);
    v=rotation.'*vParent; w=rotation.'*wParent;
    a=rotation.'*aParent; alpha=rotation.'*alphaParent; g=rotation.'*g;
end
wrench=softarm_recursive_tip(v,w,a,alpha,g,p);
tauArm=zeros(15,1);
for section=5:-1:1
    rotation=Rend(:,:,section); force=rotation*wrench(1:3); moment=rotation*wrench(4:6);
    localFirst=(section-1)*3+1; localLast=localFirst+3-1;
    tauArm(localFirst:localLast)=ownTau(:,section)+Jv(:,:,section).'*force+Jw(:,:,section).'*moment;
    wrench=ownWrench(:,section)+[force;cross(rend(:,section),force)+moment];
end
tau=[baseTau+mountMap*wrench;tauArm]+softarm_recursive_internal_force(q,dq,p);
end
