plot3dframe([0;0;0],eye(3),'k')
hold on
Gmean=mean_rigid(allG)
for(n=1:size(allG,3))
   plot3dframe(allG(1:3,4,n),allG(1:3,1:3,n),'r')
end
plot3dframe(Gmean(1:3,4),Gmean(1:3,1:3),'b')
plot3dframe(gtruth(1:3,4),gtruth(1:3,1:3),'g')
hold off
axis equal
