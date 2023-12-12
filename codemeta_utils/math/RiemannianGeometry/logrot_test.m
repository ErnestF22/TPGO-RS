Nit=10000;
werr=0;
for(it=1:Nit)
    w=randn(3,1);
    w=w/norm(w)*mod(norm(w),pi);    %make sure rotation is less than 180 degrees
    err(it)=norm(w-logrot(rot(w)));
    if(err(it)>1e-6)
        warning('Check the code!');
        return
    end
    if(err(it)>werr)
        werr=err(it);
        ww=w;
    end
end
plot(err)
norm(ww)