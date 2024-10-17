function h1=sphere_tangentProj(y,h)
flagPermute=false;
if(size(h,2)==1)
    h=permute(h,[1 3 2]);
    y=permute(y,[1 3 2]);
    flagPermute=true;
end

Nh=size(h,2);
Ny=size(y,2);

if Nh>1
    h1=zeros(size(h));
    for iN=1:Nh
        if Ny==1
            h1(:,iN)=sphere_tangentProj(y,h(:,iN));
        else
            h1(:,iN)=sphere_tangentProj(y(:,iN),h(:,iN));
        end
    end
else
    h1=h-y*y'*h;
end

if(flagPermute)
    h1=permute(h1,[1 3 2]);
end

