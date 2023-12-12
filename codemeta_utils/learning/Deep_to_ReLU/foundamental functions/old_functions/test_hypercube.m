% this funcion test if the ReluPlexFindAllVerticesOfNetwork function can 
%find all the vertices
clear all
close all
% 
d = 5;
%A=(nCons,d+1)
%A(:,1:d)x<=A(:,d+1)
% A = [ 0 1 0 7;
%     0 0 -1 -3;
%     0 -1 0 -3
%     0 0 1 7;
%     1 0 0 7;
%     -1 0 0 -3;
%     0 -2 1 0;
%     0 2 -1 9];

% A = [1 0 0 3;
%     0 1 0 3;
%     0 0 1 3;
%     -1 0 0 1;
%     0 -1 0 1;
%     0 0 -1 1];
% 
% A = [ 0 1 0 0 7;
%     0 0 -1 0 -3;
%     0 -1 0 0 -3
%     0 0 1 0 7;
%     1 0 0 0 7;
%     -1 0 0 0 -3;
%     0 -2 1 0 0;
%     0 2 -1 0 9;
%     0 0 0 1 7;
%     0 0 0 -1 -1];
% 
% A = [1 0 0 1;
%     0 1 0 1;
%     0 0 1 1;
%     0 0 -1 -1];
% 
% A = rand(4,4);
n = 9;
a = -5; b = 15;
A = rand([n,d+1]);
A(:,end) = b*rand([n,1]);
A_add = -eye(d);
A_add(:,end+1) = 0;
A = [A;A_add];


for i = 1:size(A,1)
    Aa(i,:) = A(i,1:d);
    ba(i,1) = A(i,d+1);
end

As = [Aa eye(size(Aa,1)) ba];
As(end+1,:) = 0;
% basic = [size(Aa,2)+1:size(As,2)-1];
basic = [size(Aa,2)+1:size(Aa,2)+n,1:d];
[basic,result,Ar,~]= dual_simplex2(As,basic,d);
vertices = [];
% if ~isempty(basic)
%     vertices = ReluPlexVertexDualSimplex(basic,Ar,result,d);
% end
%  plotHyperCube
% axis equal
% xlabel('x')
% ylabel('y')
% 
% for i=1:size(vertices,1)
%     plot3(vertices(i,1),vertices(i,2),vertices(i,3),'ro','Markersize',10,'MarkerFaceColor','r');
%     hold on
% end

c = 0;
 for k = 1:500
AT = A;
Q = orth(randn(d));
D = diag(randi([1,5],[1,d]));
P = Q*D*Q';
t = randi([10, 50], [d,1]);
% A(Px+t)<=b

AT(:,1:d) = AT(:,1:d)*P;
AT(:,d+1) = AT(:,d+1)+AT(:,1:d)*t;

Aa = [];
ba = [];
for i=1:size(AT,1)
    Aa(i,:) = AT(i,1:d);
    ba(i,1) = AT(i,d+1);
end

As = [Aa eye(size(Aa,1)) ba];
As(end+1,:) = 0;
basic = [size(Aa,2)+1:size(As,2)-1];
[basic,result,Ar,~] = dual_simplex2(As,basic,d);
verticesT = [];
if ~isempty(basic)
    verticesT = ReluPlexVertexDualSimplex(basic,Ar,result,d);
end
if k==1
    vertices_old = verticesT;
end
if size(verticesT,1)==size(vertices_old,1)
   c=c+1;
else
    aa=1000;
end

vertices_old = verticesT;

end

% figure(2)
% x = 0:0.2:10;
% y = 0:0.2:10;
% [X, Y] = meshgrid(x, y);
% for i=1:size(A,1)
%     Z = -(A(i,1)*X+A(i,2)*Y-A(i,4))/A(i,3);
%     surf(X, Y, Z);
%     hold on
% end
% axis equal
% xlabel('x')
% ylabel('y')
% 
% for i=1:size(verticesT,1)
%     plot3(verticesT(i,1),verticesT(i,2),verticesT(i,3),'ro','Markersize',8,'MarkerFaceColor','r');
%     hold on
% end
