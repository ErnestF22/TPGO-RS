%function [Rest,Test]=poseEstimation(X,x,varargin)
%Estimate the pose of an object having 3D shape X from the image points x
function [GEst,xest]=poseEstimationQuadratic(X,x,varargin)
methodAbsolutePoses='pose';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

x=homogeneous(x,3);

n=size(X,2);    %number of points in the model

%Compute interdistance and position of each new unknown in the solution
%vector
dij=zeros(n);
for i=1:n
    for j=i:n
        xd=X(:,i)-X(:,j);
        dij(i,j)=xd'*xd;
    end
end
[map,m]=powers2map(n);

%build the first linear system
m2=0;
for i=1:n
    for j=i+1:n
        m2=m2+1;
        M(m2,map(i,i))=x(:,i)'*x(:,i);
        M(m2,map(j,j))=x(:,j)'*x(:,j);
        M(m2,map(i,j))=-2*x(:,i)'*x(:,j);
        M(m2,m+1)=-dij(i,j);
    end
end

%find the kernel of M
[U,S,V]=svd(M);

N=n+1;
kerM=V(:,end-n:end);

%storage for optimizing repeated calls to vecsym
idxPowers=[];   %this will be filled with the first call to vecsym

%build the linear system for the lambdas
m2Lambda=0;
for i=1:n
%    for(l=1:n)
     l=i;
        for j=i+1:n
            for k=1:n
                m2Lambda=m2Lambda+1;
                vij=kerM(map(i,l),:)';
                vkl=kerM(map(j,k),:)';
                vi1j1=kerM(map(i,j),:)';
                vk1l1=kerM(map(l,k),:)';
                A=vij*vkl'-vi1j1*vk1l1';
                [K(m2Lambda,:),idxPowers]=vecsym(A,idxPowers);
            end
        end
%    end
end

%find the kernel of K (that is, veronese(lambda,2)
[Uk,Sk,Vk]=svd(K);

lambdabar=Vk(:,end);
if(lambdabar(1)<0)      %make sure x_1^2 is positive
    lambdabar=-lambdabar;
end

%find lambda
lambda=inverseVeronese(lambdabar,2,N);

%find depths
xbar=kerM*lambda;
xbar=xbar/xbar(end);

t=inverseVeronese(xbar(1:end-1),2,n);
% if(t(1)<0)
%     t=-t;
% end

t=abs(t);

%find x
xest=x.*repmat(t',3,1);

%solve absolute orientation problem
[REst,TEst]=absolute_orientation(X,xest);
GEst=RT2G(REst,TEst);

if strcmpi(methodAbsolutePoses,'reference')
    GEst=invg(GEst);
end
