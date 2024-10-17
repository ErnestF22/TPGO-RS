Nit=100;
werr=0;
for(it=1:Nit)
    u1=randn(3,1);
    u1=u1/norm(u1)*mod(norm(u1),pi/2);    %make sure rotation is less than 180 degrees
    u2=randn(3,1);
    u2=u2/norm(u2)*mod(norm(u2),pi);    %make sure rotation is less than 180 degrees
    w1=BCHrot(u1,u2);
    w2=logrot(rot(u1)*rot(u2));
    R1=rot(w1);
    R2=rot(w2);
    err(it)=acos(min(1,(trace(R1'*R2)-1)/2))*180/pi;
    if(err(it)>1e-4)
        warning('Check the code!');
        return
    end
    if(err(it)>werr)
        werr=err(it);
        uu1=u1;
        uu2=u2;
    end
end
plot(err)
%norm(ww)
