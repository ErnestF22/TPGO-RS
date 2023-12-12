function m=rot3r3_metric(g,v1,v2)
m=rot_metric(g(1:3,1:3),v1(1:3,1:3,:),v2(1:3,1:3,:));
m=m+sum(squeeze(v2(1:3,4,:).*v1(1:3,4,:)));
