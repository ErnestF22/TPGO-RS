close all
clear 
load set_out.mat
load set_all.mat

t1 = -30;
t2 = 45;

A1 = [cosd(t1) -sind(t1);sind(t1) cosd(t1)];
b1 = [1;1];

A2 = [cosd(t2) -sind(t2);sind(t2) cosd(t2)];
b2 = [1;1];

A3 = [1 -1;1 1];
b3 = [1;1];

A4 = [1 1];
b4 = 1;

%%
x=[10; -5];
xmin=x;
openfig('fig_find_minimum');
A_cum=[1.22474487139159 0.707106781186548];
b_cum=6.41421356237310;
y=A_cum*x+b_cum;
ymin=y;
figure(1);
plot3(x(1),x(2),y,'+r');

%% first
z=[1,0,1,1,1,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
% [basicd,resultd,Ard]=dual_simplex2(As,basic);
% [Ass1,basics] = getsimplexA(z,Ab_set,set_all);
% [basics,results,Ars]=simplex(Ass,basics);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');
[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);
for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
    
end
figure(1);
plot3(xmin(1),xmin(2),ymin,'+b');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[1,0,1,1,1,1,1]; 1,2
%% second [0,0,1,1,1,1,1] [0,1,1,1,1,1,1] [1,1,1,1,1,1,1]
z=[0,0,1,1,1,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');

[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);

for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
end
figure(1);
plot3(x(1),x(2),y,'+k');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[0,0,1,1,1,1,1];
%% second [0,0,1,1,1,1,1] [0,1,1,1,1,1,1] [1,1,1,1,1,1,1]
z=[0,1,1,1,1,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');

[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);

for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
end
figure(1);
plot3(x(1),x(2),y,'+k');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[0,1,1,1,1,1,1];
%% second [0,0,1,1,1,1,1] [0,1,1,1,1,1,1] [1,1,1,1,1,1,1]
z=[1,1,1,1,1,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');

[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);

for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
end
figure(1);
plot3(x(1),x(2),y,'+k');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[1,1,1,1,1,1,1]; 1,5
%% third [1,1,1,1,0,1,1] [0,1,1,1,0,1,1]
z=[1,1,1,1,0,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');

[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);

for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
end
figure(1);
plot3(x(1),x(2),y,'+k');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[1,1,1,1,0,1,1];
%% third [1,1,0,1,0,1,1] [0,1,1,1,0,1,1]
z=[0,1,1,1,0,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');

[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);

for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
end
figure(1);
plot3(x(1),x(2),y,'+k');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[0,1,1,1,0,1,1]; 1,3
%% fourth [1,1,0,1,0,1,1] [0,1,0,1,0,1,1]
z=[1,1,0,1,0,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');

[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);

for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
end
figure(1);
plot3(x(1),x(2),y,'+k');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[1,1,0,1,0,1,1];

%% fourth [1,1,0,1,0,1,1] [0,1,0,1,0,1,1]
z=[0,1,0,1,0,1,1];
[As,basic,d] = getdualAnew(z,Ab_set,set_all);
[basic,result,Ar]=dual_simplex_new(As,basic,d);
v = result(1:d)-result((d+1):2*d);
vertices=[v];
[V]=find_vertices_new(basic,Ar,d);
vertices = [vertices;V];
vertices = unique(round(vertices,4),'rows');

[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_cum=set_all.b_cum(loc,:);

for i=1:size(vertices,1)
    x=vertices(i,:)';
    y=A_cum*x+b_cum;
    if y<ymin
        disp('find new minimum')
        ymin=y
        xmin=x;
    end
end
figure(1);
plot3(x(1),x(2),y,'+g');
Az=As(1:(end-1),1:d);
bz=-As(1:(end-1),end);
yz=Az*xmin+bz; %z=[0,1,0,1,0,1,1];
