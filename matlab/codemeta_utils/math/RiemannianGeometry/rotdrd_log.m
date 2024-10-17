function v=rotdrd_log(g1,g2)
d=size(g1,2)-1;
Ng=size(g2,3);
v=zeros(size(g2));
v(1:d,1:d,:)=rot_log(g1(1:d,1:d),g2(1:d,1:d,:));
v(1:d,d+1,:)=squeeze(g2(1:d,d+1,:)-repmat(g1(1:d,d+1),[1 1 Ng]));
