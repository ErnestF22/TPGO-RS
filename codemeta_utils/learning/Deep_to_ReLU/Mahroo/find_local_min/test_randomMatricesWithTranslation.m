% this funcion test if the ReluPlexFindAllVerticesOfNetwork function can 
%find all the vertices
clear all
close all

nConstraints = 20;
rMin = 1;
rMax = 20;
d = 5;
A =  randi([rMin, rMax], [nConstraints,d+1]);

Q =  randi([-2, 2], [d,d]);
D = diag(randi([1,10],[1,d]));
P = Q*D*Q';
t = randi([rMin, rMax], [d,1]);
    
for i=1:size(A,1)
    Aa(i,:) = A(i,1:d);
    ba(i,1) = A(i,d+1);
end

As = [Aa -Aa eye(size(Aa,1)) ba; zeros(1,2*size(Aa,2)+size(Aa,1)+1)];
basic = [2*size(Aa,2)+1:size(As,2)-1];
[basic,result,Ar]= DualSimplex(As,basic,d);
vertices = [];
if ~isempty(basic)
    vertices = ReluPlexVertexDualSimplex(basic,Ar,result,d);
end

% A(Px+t)<=b
% P = diag([1:d]);
A(:,1:d) = A(:,1:d)*P;
A(:,d+1) = A(:,d+1)-A(:,1:d)*t;
Aa = [];
ba = [];
for i=1:size(A,1)
    Aa(i,:) = A(i,1:d);
    ba(i,1) = A(i,d+1);
end

As = [Aa -Aa eye(size(Aa,1)) ba; zeros(1,2*size(Aa,2)+size(Aa,1)+1)];
basic = [2*size(Aa,2)+1:size(As,2)-1];
[basic,result,Ar]= DualSimplex(As,basic,d);
verticesT = [];
if ~isempty(basic)
    verticesT = ReluPlexVertexDualSimplex(basic,Ar,result,d);
end

size(vertices,1)
size(verticesT,1)
