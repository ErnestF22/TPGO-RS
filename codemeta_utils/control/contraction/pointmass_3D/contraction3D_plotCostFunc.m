function [ ] = contraction3D_plotCostFunc( numPoints )
% Using the L1-L2 function p(x) = 2*(sqrt(1+norm(x)^2/2)-1), plot the
% eigenvalues of the hessian

%Generate random points and find the eigenvalues
eValList = [];
for i = 1:numPoints
    %Generate a random point
    x = rand(3,1)*100-50;
    %Find Hessian
    C1 = sqrt(1+norm(x)^2/2);
    H = (C1*eye(3)-1/2*x*x'/C1)/C1^2;
    %find Eigenvalues
    eG = eig(H);
    eValList = [eValList; eG];
end

%plot eigenvalues
figure
title('Eigenvalues of Hessian');
plot(eValList,'rx');

fprintf('MIN: %f, MAX: %f\n',min(eValList), max(eValList));
end

