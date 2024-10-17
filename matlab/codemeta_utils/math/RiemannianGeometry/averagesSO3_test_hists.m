%Compare different ways for averaging rotations
function averagesSO3_test_hists
clc

N=20;       % number of nodes
A=adjgallery(N,'kneigh',3);
%A=adjgallery(N,'sparserand',0.4);
%A=adjgallery(N,'full');

gtruth=exp_se3([pi;0;0],[5;5;5]); % true value on which to reach consensus
Rtruth=gtruth(1:3,1:3);
sigmanoise=0.35;                 % variance of added noise

Nit=70;    %consensus iterations
NitMode5=4; NitScalar=140;
NitMode15=2; NitScalar15=140;

Ntrials=1000;

randn('state',0)
rand('state',0)

for(itrial=1:Ntrials)
    fprintf('Trial %d/%d\n', itrial, Ntrials)
    %initialization
    for(inode=1:N)
        t_node(inode).gi=noiserigid(gtruth,sigmanoise);
        t_node(inode).aij=A(inode,:);
        t_node(inode).d=sum(t_node(inode).aij); %degree of the node
    end

    allG=get_all_rotations(t_node,'gi');
    allR=allG(1:3,1:3,:);
    
    
    err(itrial,1)=rot_dist(Rtruth,mean_rot(allR));
    err(itrial,2)=rot_dist(Rtruth,mean_rot_quat(allR));
    err(itrial,3)=rot_dist(Rtruth,mean_rot_eucl(allR));
end

legendtext=char('Riemannian','Quaternions','Euclidean');

save averagesSO3_test_hists_results

function R=mean_rot_quat(allR)
q=rot2quat(allR);
q=quat_align(q);
R=quat2rot(sum(q,2));

function R=mean_rot_eucl(allR)
R=rot_proj(sum(allR,3));
