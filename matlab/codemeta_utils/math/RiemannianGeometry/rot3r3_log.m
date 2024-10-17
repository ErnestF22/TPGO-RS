function v=rot3r3_log(g1,g2)
Ng=size(g2,3);
v=zeros(size(g2));
v(1:3,1:3,:)=rot_log(g1(1:3,1:3),g2(1:3,1:3,:));
v(1:3,4,:)=squeeze(g2(1:3,4,:)-repmat(g1(1:3,4),[1 1 Ng]));
