function hessDistance_test
wi=[1 1 0];
wj=[1 0 1];

wi=wi/norm(wi);
wj=wj/norm(wj);

Np=100;
x=linspace(-pi,pi,Np);

[X,Y]=meshgrid(x,x);
Z=zeros(size(X));

for(ix=1:Np)
    for(iy=1:Np)
        Z(ix,iy)=hessDistance(X(ix,iy)*wi,Y(ix,iy)*wj);
    end
end

figure(1)
surf(X,Y,Z,'EdgeColor','none')
set(gcf,'Renderer','zbuffer')
