function [V,T,R] = solving_dual_simplex(R)
A=[];
A = [-R(:,1:2) ones(size(R,1),1) eye(size(R,1)) R(:,3)];
A = [A; zeros(1,size(A,2))];

basic = 6:5+size(R,1);
[basic,~,A,T]=dual_simplex(A,basic);
V = result(A,basic);
[V,R] = find_vertices(A,basic);

for i=1:size(V,1)
    plot3(V(i,1),V(i,2),V(i,3),'r*')
    hold on
end

end