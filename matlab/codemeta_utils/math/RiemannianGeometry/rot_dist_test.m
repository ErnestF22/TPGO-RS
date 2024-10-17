function rot_dist_test
Sij=rot([1;0;0]);
wiwj1=[-1 -1 1 -1 -1 0];
wiwj2=[-1 -1 0 -1 -1 1];

wiwj1=wiwj1/norm(wiwj1);
wiwj2=wiwj2/norm(wiwj2);

Np=100;
x=linspace(-pi,pi,Np);

[X,Y]=meshgrid(x,x);
Z=zeros(size(X));

for(ix=1:Np)
    for(iy=1:Np)
        wi=X(ix,iy)*wiwj1(1:3)+Y(ix,iy)*wiwj2(1:3);
        wj=X(ix,iy)*wiwj1(4:6)+Y(ix,iy)*wiwj2(4:6);
        
        Z(ix,iy)=rot_dist(rot(wi),Sij*rot(wj))^2;
    end
end

figure(1)
surf(X,Y,Z,'EdgeColor','none')
set(gcf,'Renderer','zbuffer')
