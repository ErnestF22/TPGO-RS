close all
clear
addpath('danyang/')

t1 = -30;
t2 = 45;

A1 = [cosd(t1) -sind(t1);sind(t1) cosd(t1)];
b1 = [1;1];

A2 = [cosd(t2) -sind(t2);sind(t2) cosd(t2)];
b2 = [1;1];

A3 = [1 1];
b3 = 1;

nbColors=32;
colors=parula(nbColors);

%% 
% number of layers
% L = 4;
L = 3;
% generate all possible z
set_all = table;
n = 0;
for i=1:L
    A = eval(['A' num2str(i)]);
    n = n+size(A,1);
end
z_set=dec2bin(0:2^n-1)-'0';
set_all.z = z_set;
% generate A_out and b_out for z, y=A_out*x+b_out
Ab_set = struct; % Ab_set store A's and b's at each layer
for i=1:L
    Ab_set(i).A = eval(['A' num2str(i)]);
    Ab_set(i).b = eval(['b' num2str(i)]);
end
m = size(z_set,1);
A_out=[];
b_out=[];
for i=1:m
    z = z_set(i,:)';
    [A,b] = get_cumulative(Ab_set,L,z);
    A_out = [A_out; A];
    b_out = [b_out; b];
end
set_all.A_out = A_out;
set_all.b_out = b_out;

%% exclude 0 rows
Ab = [A_out,b_out];
idx = find(sum(Ab,2)==0);
set_all(idx,:) = [] ;

%% find the neighborhood region
z_set = set_all.z;
legend_z=[];
set_out = struct;
j=1;
% x=linspace(-15,15,100);
x = [1;1];
% zset=[1,1,1,1,1; 1,1,0,1,1; 1,0,1,1,1; 0,1,0,1,1; 0,1,1,1,1; 0,0,1,1,1];
zset=[1,1,1,1,1; 1,1,0,1,1]';
Q1 = [];
for z = zset
    [As,basic] = getdualAmatrix(z,Ab_set);
    [basic,Ar,P]=dual_simplex_xpn(As,basic);
    if isempty(basic)
        continue
    else
        [Q,ACT] = find_vertices(basic,Ar,z',x);
        if isempty(Q1)
            Q1=Q;
        else
            Q1 = [Q1,Q];
        end
    end
end
test = 1;
%% check flip matrix
% result = check_tableau_after_flipping([1,0,1,1,1]',1,P,As,Ar,Ab_set,active_idx)
%% find the neighborhood region
z_set = set_all.z;
legend_z=[];
set_out = struct;
j=1;
x=linspace(-15,15,100);
% for i=1:size(z_set,1)
%     z = z_set(i,:)';
for i=1:1
    z=[1,1,1,1,1]';
%     z=[1,1,0,1,1]';
%     z=[1,0,1,1,1]';
%     z=[0,1,0,1,1]';
%     z=[0,1,1,1,1]';
%     z=[0,0,1,1,1]';
    [As,basic] = getdualAmatrix(z,Ab_set);
    [basic,Ar,P]=dual_simplex_xpn(As,basic);
    d = (size(Ar,2)-size(Ar,1))/2;
    if isempty(basic)
        continue
    else
        % if the region exist
        set_out(j).z=z';
        % transform to Az*x+bz<=0, Az and bz are cumulative A and b of each
        % layer, y1=Az(1,:)*x+bz(1,:)
        set_out(j).Az=As(1:(end-1),1:d);
        set_out(j).bz=-As(1:(end-1),end);
        j = j+1;
        v = result(1:d)-result((d+1):2*d);
        [ACT, V]=find_vertices_query_xpn(basic,Ar,2*d);
        visited_V = find_vertices_all_points(basic,Ar,2*d);
        vertices = [v;V];
    end
end
vertices = unique(round(vertices,4),'rows')

%% plot All A's and b's
ax = [-15 15;-15 15];
A_out=[];
b_out=[];
for i=1:size(set_out,2)
    Az = set_out(i).Az;
    bz = set_out(i).bz;
    for j=1:size(Az,1)
        A_out=[A_out;Az(j,:)];
        b_out=[b_out;bz(j,:)];
    end
end
figure(1);
plot2doutput(A_out,b_out,ax)
hold on

%% plot vertices
figure(1);
for i=1:size(vertices,1)
    plotPoints(vertices(i,:)','r','Markersize',10);
    hold on;
end

%% plot visited points
figure(1);
for i=1:size(vertices,1)
    plotPoints(visited_V(i,:)','k','Markersize',10);
    hold on;
end

%% construct whole output space
x=linspace(-15,15,100);
ax = [-15 15;-15 15];
for j=1:size(set_out,2)
    A_out=[];
    b_out=[];
    z=set_out(j).z;
    for i=1:size(vertices,1)
        vertex = vertices(i,:);
        idx = find_z(set_out,z);
        Az = set_out(idx).Az;
        bz = set_out(idx).bz;
        out = round(Az*vertex'+bz,4);
        if sum(out==0)>=2
            idx=find(out==0);
            A_out=[A_out;Az(idx,:)];
            b_out=[b_out;bz(idx,:)];
        end
    end
    [~,loc] = ismember(z,set_all.z,'rows');
    A_cum=set_all.A_cum(loc,:);
    b_cum=set_all.b_cum(loc,:);
    figure(1);
    funImage(x,x,@(x) halfspaceout(x,A_out,b_out,A_cum,b_cum),...
        'method','surf','methodOpts',{'FaceColor',colors(mod(j*n,nbColors)+1,:)});
    hold on
end
% legend('{[1,1,1,1,1]}','{[1,1,0,1,1]}','{[1,0,1,1,1]}','{[0,1,1,1,1]}','{[0,1,0,1,1]}','{[0,0,1,1,1]}')
% legend('1,1,1,1,1','1,1,0,1,1','1,0,1,1,1','0,1,1,1,1');
% legend('0,1,1,1,1','1,1,1,1,1','0,0,1,1,1','0,1,0,1,1');
% legendCell = cellstr(num2str(legend_z));
% legend(legendCell)

%% construct output space for one z
z=[1,1,1,1,1];
% z=[1,1,0,1,1];
% z=[1,1,1,1,0,1,1];
% z=[1,0,1,1,0,1,1];
ax = [-15 15;-15 15];
x=linspace(-15,15,100);
A_out=[];
b_out=[];
for i=1:size(vertices,1)
    vertex = vertices(i,:);
    idx = find_z(set_out,z);
    Az = set_out(idx).Az;
    bz = set_out(idx).bz;
    out = round(Az*vertex'+bz,4);
    if sum(out==0)>=d
        idx=find(out==0);
        A_out=[A_out;Az(idx,:)];
        b_out=[b_out;bz(idx,:)];
    end
%     A_out=[A_out;Az];
%     b_out=[b_out;bz];
end
[~,loc] = ismember(z,set_all.z,'rows');
A_cum=set_all.A_cum(loc,:);
b_out=set_all.b_cum(loc,:);
figure(1);
x=linspace(-15,15,100);
funImage(x,x,@(x) halfspaceout(x,A_out,b_out,A_cum,b_out),...
    'method','surf','methodOpts',{'FaceColor',colors(mod(randi(32),nbColors)+1,:)});
hold on
% plot2doutput(A_out,b_out,ax);

%% choose a vertex
vertex = vertices(1,:);
z=[];
for i=1:size(set_out,2)
    Az = set_out(i).Az;
    bz = set_out(i).bz;
    zi = set_out(i).z';
    out = round(Az*vertex'+bz,4);
    if sum(out==0)>=2
        z=[z,zi'];
    end
end

%% construct All A's and b's
z=[];
ax = [-15 15;-15 15];
A_out=[];
b_out=[];
for i=1:size(vertices,1)
    vertex = vertices(i,:);
    for j=1:size(set_out,2)
        Az = set_out(j).Az;
        bz = set_out(j).bz;
        zi = set_out(j).z';
        out = round(Az*vertex'+bz,4);
        if sum(out==0)>=2
            idx=find(out==0);
            A_out=[A_out;Az(idx,:)];
            b_out=[b_out;bz(idx,:)];
        end
    end
end
figure(1);
plot2doutput(A_out,b_out,ax)
hold on

%%
x=[10; -5];
openfig('fig_find_minimum');
A_cum=[1.22474487139159 0.707106781186548];
b_out=6.41421356237310;
y=A_cum*x+b_out;
figure(1);
plot3(x(1),x(2),y,'+r');

%%
x=[10; -5];
openfig('layer3fig');