close all
clear
addpath('foundamental functions/')
addpath('find_plane_v2')

t1 = -30;
t2 = 45;

A1 = [cosd(t1) -sind(t1);
    sind(t1) cosd(t1);
    cosd(t2) -sind(t2);];
b1 = [1;1;1];

A2 = [cosd(t2) -sind(t2) 0;
    sind(t2) cosd(t2) 0];
b2 = [1;1];

A3 = [-10 1];
b3 = 0.5;

nbColors=128;
colors=parula(nbColors);

% number of layers
L = 3;

% generate all possible z
n = 0;
for i=1:L
    A = eval(['A' num2str(i)]);
    n = n+size(A,1);
end
z_set=dec2bin(0:2^n-1)-'0';

for i=1:L
    Ab_set(i).A = eval(['A' num2str(i)]);
    Ab_set(i).b = eval(['b' num2str(i)]);
end

%% plot output space
figure(1);
x = -10:0.1:10;
z = {};
nbColors=64;
colors=parula(nbColors);
gs1=getSeqCombSet([0 1],3);
cnt=0;
for is1=1:8
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
legend(legendInfo);
savefig('outputspace.fig');
%%
close all;
openfig('outputspace.fig');
hold on;
% plot input point
x = [10;-10];
plotPoints(x,'k','Markersize',10);
hold on;

ax = [-10 10;-10 10];
[~,z] = forward(x,L,Ab_set); % get z with input x
[As,basic] = getdualAmatrix(z,Ab_set);
d = size(As,2)-size(As,1);
[basic,Ar,P] = dual_simplex_xpn(As,basic);
if ~isempty(basic) % if basic is empty -> no violating point
    [Q,V,ZSet] = find_plane_P(basic,As,Ar,P,z,x,Ab_set);
else
    disp('No violating point found');
    Q = [];
    min_violate_idx = [];
end
% plot phase of first region
for j = 1:size(Q,2)
    planeidx = Q(j).planeidx-d;
    Al = As(planeidx,1:d/2);
    bl = -As(planeidx,end);
    plot2doutput(Al,bl,ax);
    hold on;
end

i = 1;
flagCheck = false;
d = size(As,2)-size(As,1);
while ~isempty(Q) % if there are unchecked points
    h_plane_set = extractfield(Q,'h_plane'); % get the plane distance list
    [min_distance,min_idx] = min(h_plane_set); % find minimum distance
    
    z = Q(min_idx).z;
    zflip = Q(min_idx).planeidx-d;
    z(zflip) = ~z(zflip);
    if (sum(sum(abs(ZSet-z'),2)==0)==0) % if z not checked
        [Q,ZSet,act_set,As,V] = find_neighbor_plane(Q,min_idx,x,ZSet,d,Ab_set);
        close all;
        openfig('outputspace.fig');
        hold on;
        for k = 1:size(V,1)
            vertex = V(k,:);
            plotPoints(vertex','r','Markersize',10);
            hold on;
        end
        for j = 0:(size(act_set,2)-1)
            planeidx = Q(end-j).planeidx-d;
            Al = As(planeidx,1:d/2);
            bl = -As(planeidx,end);
            plot2doutput(Al,bl,ax);
            hold on;
        end
    end
    C(i) = Q(min_idx); % put it to C
    Q(min_idx)=[]; % pop out from Q
    i = i+1;
    flagCheck = false;
    if mod(i,1)==0
        checkQueue = ['round ',num2str(i),', waiting for checking: ',num2str(size(Q,2))];
        disp(checkQueue);
        distancei = ['min distance is ',num2str(min_distance)];
        disp(distancei);
    end
end