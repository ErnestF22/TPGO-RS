%Finds a and x such that kerM*a=veronese(x,2)
%function [x,a]=solveVeronese2LinearCombination(kerM,dimX)
%Compute the intersection of the range of kerM and the image of veronese(x,2)
%Inputs
%   kerM    matrix whose range contains the veronese map of x
%   dimX    dimension of x
function [x,lambda,K,lambdabar]=solveVeronese2LinearCombination(kerM,dimX)

n=dimX;
map=powers2map(n);

%storage for optimizing repeated calls to vecsym
idxPowers=[];   %this will be filled with the first call to vecsym
K=[];

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

%check that range(kerM) is not too large
if rank(K)<(size(K,2)-1)
    error('The range of kerM is too large to find a unique solution (i.e., kerM has too many columns')
end

%find the kernel of K (that is, veronese(lambda,2)
[Uk,Sk,Vk]=svd(K);

lambdabar=Vk(:,end);

%find a
if(lambdabar(1)<0)      %make sure x_1^2 is positive
    lambdabar=-lambdabar;
end
lambda=inverseVeronese(lambdabar,2,size(kerM,2));

%find x
xbar=kerM*lambda;
if(xbar(1)<0)      %make sure x_1^2 is positive
    xbar=-xbar;
end
x=inverseVeronese(xbar,2,dimX);
