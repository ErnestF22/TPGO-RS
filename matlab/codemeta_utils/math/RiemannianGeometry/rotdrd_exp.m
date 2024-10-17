function g=rotdrd_exp(g,v)
d=size(g,2)-1;
g(1:d,1:d,:)=rot_exp(g(1:d,1:d,:),v(1:d,1:d,:));
g(1:d,d+1,:)=g(1:d,d+1,:)+v(1:d,d+1,:);
