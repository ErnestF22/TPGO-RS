function [X,Xn,E] =  CreateData(sigma) 




%% load data

load joints15
load('15_02.mat')

Sl = S(:,mocapIndices);


%% create edge set
E = [];
for i=1:numel(skel.tree)
    for j=1:numel(skel.tree(i).children)
       E = [E;i skel.tree(i).children(j) ] ;
    end  
end

[~,idx] = sort(E(:,2));
E = E(idx,:);

%% create noisy data
F = size(Sl,1)/3; % # of frames
P = size(Sl,2); % # of points/joints




X = zeros([3 F P]); % groundtruth
Xn = zeros(size(X)); % noisy
for p=1:P, 
    X(:,:,p) = reshape(Sl(:,p),[3 F]); 
    Xn(:,:,p) = X(:,:,p) + sigma*randn(size(X(:,:,p)));
end

X = permute(X,[1 3 2]); % 3 x P x F
Xn = permute(Xn,[1 3 2]); % 3 x P x F

end
