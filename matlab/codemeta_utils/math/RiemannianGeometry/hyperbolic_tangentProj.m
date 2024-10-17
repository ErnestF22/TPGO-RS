function v=hyperbolic_tangentProj(x,v)
flagPermute=false;
if(size(v,2)==1)
    v=permute(v,[1 3 2]);
    flagPermute=true;
end

xp=cnormalize([x(1:end-1);-x(end)]);

[D,N]=size(v);

ip=sum(xp(:,ones(1,N)).*v);

v=v-ip(ones(1,D),:).*xp(:,ones(1,N));

if(flagPermute)
    v=permute(v,[1 3 2]);
end
