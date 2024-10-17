function POCVisualizationSphereLog
k=5;
[x,y,z] = sphere(2^k-1);
[xSub,ySub,zSub] = sphere(16);
y0=cnormalize([0;-1;1]);
f=@(x) sphere_dist(x,y0);
gradf=@(x) -sphere_log(x,y0);
xyz=cat(3,x,y,z);
xyzSub=cat(3,xSub,ySub,zSub);
c=evalfunVec(f,xyz);
g=-evalfunVec(gradf,xyzSub);

figure(1)
surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),c,'EdgeColor','Interp','FaceColor','Interp');
hold on
quiver3(xyzSub(:,:,1),xyzSub(:,:,2),xyzSub(:,:,3),g(:,:,1),g(:,:,2),g(:,:,3),1.5,'k','LineWidth',2)
a=0.5;
colormap(jet*a+(1-a))
hold off
axis equal
