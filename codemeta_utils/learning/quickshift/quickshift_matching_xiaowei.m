function [membership,info] = quickshift_matching(D,groupLabel,varargin)

rho = 1;
threshold=0.01;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'rho'
            ivarargin=ivarargin+1;
            rho=varargin{ivarargin};
        case 'threshold'
            ivarargin=ivarargin+1;
            threshold=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%minimal distance
w = ones(size(D,1),1);
for i = unique(groupLabel)
    ind = find(groupLabel==i);
    B = D(ind,ind);
    B(1:size(B,1)+1:end) = inf;
    w(ind) = min(B);
end

%density
P = sum(bsxfun(@times,w,exp(-bsxfun(@rdivide,D.^2,rho*(w+eps).^2))),1);

[vDist,vMember]=quickshift_tree(P,D);
vMemberTh=quickshift_breakTree(vDist,vMember,'threshold',threshold);
membership=quickshift_tree2membership(vMemberTh);

if nargout>1
    info.density=P;
    info.treeVectorMembership=vMember;
    info.treeVectorMembershipThresholded=vMemberTh;
    info.treeVectorDistances=vDist;
end
