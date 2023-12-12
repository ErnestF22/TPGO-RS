% this funcion test if the ReluPlexFindAllVerticesOfNetwork function can 
%find all the vertices
clear all
close all
addpath('old_functions')
% A(:,1:d)x<=A(:,d+1)
% all vertices(generate randomly) in any space
% x+,x-

d = 3;

n = 10;
a = -5; b = 15;
A = a + (b-a) * rand([n,d]);
A(:,end+1) = b*rand([n,1]);
% 
% A = rand([n,d]);
% A(:,end+1) = rand([n,1]);

% load counterexample case(if stored)
% load('At.mat');
% load('Atset.mat');
% load('A_v21.mat');
% load('FailMat1.mat');
A = A;
for i = 1:size(A,1)
    Aa(i,:) = A(i,1:d);
    ba(i,1) = A(i,d+1);
end

As = [Aa -Aa eye(size(Aa,1)) ba];
As(end+1,:) = 0;
basic = [2*d+1:size(As,2)-1];
% [basic,result,Ar,~]= dual_simplex_xpn(As,basic,2*d);
[basic,result,Ar,~]= dual_simplex_xpn_old(As,basic,2*d);

vertices = [];
acts = [];
if ~isempty(basic)
    idx = (2*d+1:size(Ar,2)-1);
    a = idx(~ismember(idx,basic));
    a = sort(a,'ascend');
    v = result(1:d)-result(d+1:2*d);
    [ACT, V] = find_vertices_query_xpn_old(basic,Ar,2*d);
%     [ACT1, V1] = find_vertices_query_xpn(basic,Ar,2*d);
    vertices = [V;v];
    acts = [ACT;a];
    vertices = unique(round(vertices,4),'rows');
    acts = unique(acts,'rows');
end

c = 0;
UnsuccessfulCase = [];
Atset = {}; % store counterexample case
 for k = 1:500
AT = A;
Q = orth(randn(d));
D = diag(randi([1,5],[1,d]));
P = Q*D*Q';
t = randi([10,50], [d,1]);
% A(Px+t)<=b
AT(:,1:d) = AT(:,1:d)*P;
AT(:,d+1) = AT(:,d+1)-AT(:,1:d)*t;
Atset{k} = AT; % store counterexample case
% AT = Atset{k}; % load counterexample case
Aa = [];
ba = [];
for i=1:size(AT,1)
    Aa(i,:) = AT(i,1:d);
    ba(i,1) = AT(i,d+1);
end

As = [Aa -Aa eye(size(Aa,1)) ba];
As(end+1,:) = 0;
basic = [2*size(Aa,2)+1:size(As,2)-1];
[basic,A,P]= dual_simplex_xpn(As,basic);
% [basic,result,Ar,~]= dual_simplex_xpn_old(As,basic,2*d);
verticesT = [];
actT = [];
if ~isempty(basic)
    idx = (2*d+1:size(Ar,2)-1);
    a = idx(~ismember(idx,basic));
    a = sort(a,'ascend');
    v = result(1:d)-result(d+1:2*d);
    [ACT, V] = find_vertices_query_xpn(basic,Ar);
%     [ACT1, V1] = find_vertices_query_xpn_old(basic,Ar,2*d);
    verticesT = [V;v];
    actT = [ACT;a];
    uniqueverticesT = unique(round(verticesT,4),'rows');
    actT = unique(actT,'rows');
end
% if size(uniqueverticesT,1)==size(vertices,1)
if size(actT,1)==size(acts,1)
   c=c+1;
else
    aa=1000;
    sizeactT = size(actT,1);
    UnsuccessfulCase = [UnsuccessfulCase;k,sizeactT,isempty(basic)]; % UnsuccessfulCase = [k, number of vertices find, if the basic is empty from dual_simplex]
end
 end

% save('At.mat','A');  % store counterexample
% save('Atset.mat','Atset'); % store counterexample