close all
% clear
addpath('foundamental functions/')
addpath('nnet/')
addpath('old_version_functions')

load('/Users/danyangli/Documents/mip/Ab_set.mat');
load('/Users/danyangli/Documents/mip/x_mip.mat');
load('/Users/danyangli/Documents/mip/y_mip.mat');
load('/Users/danyangli/Documents/mip/z_mip.mat');
load('x_input.mat');
load('v_simplex.mat');
load('z_simplex.mat');
% difference between MIP and simplex
z_mip = z_mip';
% find(z_mip~=z_simplex)

%% check plane from MIP
[As,~] = getdualAmatrix(z_mip,Ab_set);
d = size(As,2)-size(As,1);
A1 = As(1:(end-1),1:d/2);
b1 = As(1:(end-1),end);
result = A1*x_mip'-b1;
find(round(result,1)==0)

%% check distance from simplex
z = z_mip;
[As,basic] = getdualAmatrix(z,Ab_set);
[basic,Ar,P] = dual_simplex_xpn(As,basic);
[ACT, vertices]=find_vertices_query_xpn(basic,Ar);
plane = unique(ACT)
if ~isempty(basic) % if basic is empty -> no violating point
    [Q,PlaneSet,ZSet] = find_plane_P(basic,As,Ar,P,z,x_input,Ab_set);
end
h_plane_set = extractfield(Q,'h_plane'); % get the plane distance list
[min_distance,min_idx] = min(h_plane_set); % find minimum distance
min_distance

%% 
figure(1);
x = -10:0.1:10;
z = {};
nbColors=32;
colors=parula(nbColors);
gs1=getSeqCombSet([0 1],2);
cnt=0;
for is1=1:4
    z1=gs1()';
    zs(1).z = z1;
    gs2=getSeqCombSet([0 1],2);
    for is2=1:4
        z2=gs2()';
        zs(2).z = z2;
        gs3=getSeqCombSet([0 1],1);
        for is3=1:2
            z3=gs3()';
            zs(3).z = z3;
            cnt = cnt+1;
            z = zs(1).z;
            for i = 2:3
                z = [z;zs(i).z];
            end
            Y = funImage(x,x,@(x) netForward(Ab_set,z,x),'method','surf','methodOpts',{'FaceColor',colors(mod(cnt-1,nbColors)+1,:)});
            if ~all(all(isnan(Y)))
                legendInfo{cnt} = ['Z = ' num2str(z')];
            else
                legendInfo{cnt} = ('');
            end
            hold on
            
        end
    end
end

% plotPoints(x_input,'g','Markersize',10);
% legendInfo{end+1} = ('input point');
% hold on;
% plotPoints(x_mip','r','Markersize',10);
% legendInfo{end+1} = ('v from mip');
% hold on;
% plotPoints(v_simplex','b','Markersize',10);
% legendInfo{end+1} = ('v from simplex');

legend(legendInfo)

%% 
figure(1);

Ab_set = {};
Ab_set(1).A = [];
Ab_set(1).b = [];
L = 4;
for i =1:(L-1)
    d = 2;
    limL = -15; limU = 15;
    A1 = limL + (limU-limL) * rand([d,d]);
    b1 = limU*rand(2,1);
    Ab_set(end+1).A = A1;
    Ab_set(end).b = b1;
end
A1 = limL + (limU-limL) * rand([1,d]);
b1 = limU*rand(1,1);
Ab_set(end+1).A = A1;
Ab_set(end).b = b1;
Ab_set(1) = [];

x = -10:0.1:10;
z = {};
nbColors=129;
colors=parula(nbColors);
gs1=getSeqCombSet([0 1],2);
cnt=0;
for is1=1:4
    z1=gs1()';
    zs(1).z = z1;
    gs2=getSeqCombSet([0 1],2);
    for is2=1:4
        z2=gs2()';
        zs(2).z = z2;
        gs3=getSeqCombSet([0 1],2);
        for is3=1:4
            z3=gs3()';
            zs(3).z = z3;
            gs4=getSeqCombSet([0 1],1);
            for is4 = 1:2
                z4=gs4()';
                zs(4).z = z4;
                cnt = cnt+1;
                z = zs(1).z;
                for i = 2:4
                    z = [z;zs(i).z];
                end
                Y = funImage(x,x,@(x) netForward(Ab_set,z,x),'method','surf','methodOpts',{'FaceColor',colors(mod(cnt-1,nbColors)+1,:)});
                if ~all(all(isnan(Y)))
                    legendInfo{cnt} = ['Z = ' num2str(z')];
                else
                    legendInfo{cnt} = ('');
                end
                hold on
            end
        end
    end
end

plotPoints(x_input,'g','Markersize',10);
legendInfo{end+1} = ('input point');
hold on;
% plotPoints(x_mip','r','Markersize',10);
% legendInfo{end+1} = ('v from mip');
% hold on;
% plotPoints(v_simplex','b','Markersize',10);
% legendInfo{end+1} = ('v from simplex');

legend(legendInfo);
%%
L = size(Ab_set,2);
[y,z] = forward(x_mip',L,Ab_set); % get z with input
[As,basic] = getdualAmatrix(z,Ab_set);
[basic,Ar,P] = dual_simplex_xpn(As,basic);
plane.A = Ab_set(end).A;
plane.b = Ab_set(end).b;
if ~isempty(basic)
    [Q,ACT] = find_vertices_P(basic,Ar,P,z,x_input,plane);
else
    disp('No vertices found');
end

v_size = size(Q(1).v,2);
vset = reshape(extractfield(Q,'v'),[v_size,size(Q,2)]);
vset = vset';
[~,loc] = ismember(x_mip,vset,'rows');