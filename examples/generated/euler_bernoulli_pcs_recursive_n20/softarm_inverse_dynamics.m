function tau=softarm_inverse_dynamics(q,dq,ddq,p)
%SOFTARM_INVERSE_DYNAMICS Recursive section inverse dynamics.
%#codegen
assert(numel(q)==40&&numel(dq)==40&&numel(ddq)==40);
q=q(:); dq=dq(:); ddq=ddq(:);
base=softarm_recursive_base(q,dq,ddq(1:0),p);
assert(numel(base)==15);
v=base(1:3); w=base(4:6); a=base(7:9); alpha=base(10:12); g=base(13:15);
baseTau=base(16:15+0);
mountMap=reshape(base(16+0:end),0,6);
Rend=zeros(3,3,20); rend=zeros(3,20);
Jv=zeros(3,2,20); Jw=zeros(3,2,20);
ownWrench=zeros(6,20); ownTau=zeros(2,20);
for section=1:20
    first=0+(section-1)*2+1; last=first+2-1;
    jet=softarm_recursive_end(section,q,dq,p);
    cursor=1; Rend(:,:,section)=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    rend(:,section)=jet(cursor:cursor+2); cursor=cursor+3;
    Jv(:,:,section)=reshape(jet(cursor:cursor+5),3,2); cursor=cursor+6;
    Jw(:,:,section)=reshape(jet(cursor:cursor+5),3,2); cursor=cursor+6;
    Jvd=reshape(jet(cursor:cursor+5),3,2); cursor=cursor+6;
    Jwd=reshape(jet(cursor:cursor+5),3,2);
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
tauArm=zeros(40,1);
for section=20:-1:1
    rotation=Rend(:,:,section); force=rotation*wrench(1:3); moment=rotation*wrench(4:6);
    localFirst=(section-1)*2+1; localLast=localFirst+2-1;
    tauArm(localFirst:localLast)=ownTau(:,section)+Jv(:,:,section).'*force+Jw(:,:,section).'*moment;
    wrench=ownWrench(:,section)+[force;cross(rend(:,section),force)+moment];
end
tau=[baseTau+mountMap*wrench;tauArm]+softarm_recursive_internal_force(q,dq,p);
end
