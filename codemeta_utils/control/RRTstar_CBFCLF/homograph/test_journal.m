clear
close all

%% Environment
nodes = [15 30;20 60;45 80;90 60;70 13;30 0;30 40;40 55;50 58;65 50;60 35;40 30];
sets = [ 8 7 1 2;
         9 8 2 3;
        10 9 3 4;
        11 10 4 5;
        12 11 5 6;
        7 12 6 1];
    
landmarks(1).L = [10 50;50 50;30 80;50 40;85 40]';
landmarks(2).L = [10 50;50 50;90 70;50 0;10 20]';
for i=1:size(sets,1)
    node_set = sets(i,:);
    y = nodes(node_set,:)';
    % plot_convex(Ax,bx)
    p = (y(:,2)+y(:,4))/2;
    [A,b] = convexSet(y);
    [Ax,bx] = LPSet(A,b,p);
    Ah = -Ax;
    bh = bx;
    Ah(4,:)=[];
    bh(4)=[];
    exitDir =  -Ax(4,:)';
    exitDir = exitDir/norm(exitDir);
    % xe = intersectionTwoLines(A(1,1),b(1),A(3,1),b(3));
    xe = (y(:,1)+y(:,4))/2;
    M(i).y = y;
    M(i).Ax = Ax;
    M(i).bx = bx;
    M(i).Ah = Ah;
    M(i).bh = bh;
    M(i).xe = xe;
    M(i).xc = [y(:,1) y(:,4)];
    M(i).z = exitDir;
    %   M(i).L = y+5*rand(2,4)-10*ones(2,4);
    % if i==4
    %     M(i).L = [10 50;50 50;90 70;50 0;85 40]';
    % else
    M(i).L = landmarks(fix((i-1)/3)+1).L;%[10 50;50 50;30 80;70 80;85 40;10 50;50 50;90 40;50 0;10 20]';%;nodes';%
    % end
    
end
%% Compute control gain with new version of cost function (journal) 
% M(1).L = [5 45;35 45;25 70;45 50]';
% M(2).L = [45 50;25 70;52 50;65 75]';
% M(3).L = [52 50;65 75;60 42;85 33]';
% M(4).L = [60 42;85 33;60 2;45 53]';
% M(5).L = [60 2;45 53;18 10;36 36]';
% M(6).L = [18 10;36 36;5 45;35 45]';

start = [22;40];
flag_loop = 1;

CB = 1*[1 1 1 1 1 1 1];
CL = [1 1 1 1 1 1];
for i=1:size(sets,1)
    M(i).Cb = CB(i);
    M(i).Cl = CL(i);
end
K_new_2 = optFirstorderU_newCost_journal(M,0.1,flag_loop);
plot_journal(M,K_new_2,start,flag_loop,1)

CB = 1*[1 1 1 1 1 1 1];
CL = [2 2 2 2 1 1];
for i=1:size(sets,1)
    M(i).Cb = CB(i);
    M(i).Cl = CL(i);
end
K_new_7 = optFirstorderU_newCost_journal(M,0.7,flag_loop);
plot_journal(M,K_new_7,start,flag_loop,2)

CB = 5*[1 1 1 1 1 1 1];
CL = [1 1 1 1 1 1];
for i=1:size(sets,1)
    M(i).Cb = CB(i);
    M(i).Cl = CL(i);
end
K_new_P = optFirstorderU_newCost_journal(M,1,flag_loop);
plot_journal(M,K_new_P,start,flag_loop,3)


%% Compute control gain with previous optimization problem
Wb = 1*ones(size(Ah,1),1); 
Wl = 1;
Cb = 0.5;
Cl = 1;
for i=1:size(M,2)
    % here:
    [K,~,~,~,~,k_added] = optFirstorderWithU(Wb,Wl,Cb,Cl,M(i).z,M(i).Ah,M(i).Ax,M(i).bx,M(i).bh,M(i).L,M(i).xe,i);
    K_old(i).K = K;
    K_old(i).added = k_added;
end
plot_journal(M,K_old,start,flag_loop,4)

% navigate_journsl(start,K_old,M,'r*');
% title('red: old, blue: P=0.5, purple:P=1')
% axis equal
