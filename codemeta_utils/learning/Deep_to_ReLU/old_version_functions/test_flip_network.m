close all
clear
addpath('danyang/')
addpath('nnet/')
addpath('old_version_functions')
% this function test ACASXU network, tranform once

i = 1;
j = 3;
network = ['nnet/ACASXU_run2a_' num2str(i) '_' num2str(j) '_batch_2000.nnet'];
Ab_set = load_network(network);

x = [55947.691;0;0;1145;60];
L = size(Ab_set,2);
[~,z] = forward(x,L,Ab_set); % get z with input x
[As,basic] = getdualAmatrix(z,Ab_set);
[basic,Ar,P] = dual_simplex_xpn(As,basic);
if ~isempty(basic)
    [Q,ACT] = find_vertices_P(basic,Ar,P,z,x);
else
    Q = [];
    ACT = [];
end

% check flipped matrix
Ad = Q(1).A;
d = size(Ad,2)-size(Ad,1);
for i = 1:size(Q,2)
% for i=45
    a = Q(i).act;
    z = Q(i).z;
    P = Q(i).P;
    A = Q(i).A;
    zflip_set = unique(a)-d;
    for zflip = zflip_set
        [A1,~] = getdualAmatrix(z,Ab_set);
        [Arf,~] = find_tableau_after_flip_P(z,zflip,P,A1,A,Ab_set);
        [Q1,ACT1] = checkA(z,zflip,Ab_set,x);
        [~,col] = ismember(a,ACT1,'rows'); % find the position of z in Q
        sumall = abs(sum((Q1(col).A-Arf),'all'));
        if sumall<=1e-3
            dis = ['round ',num2str(i),' match!'];
            disp(dis);
        else
            dis = ['round ',num2str(i),' not match!'];
            disp(dis);
            disp(sumall);
        end
    end
end



