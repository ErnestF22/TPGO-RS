function [h] = contraction3D_plotEval(eMin, eMax, kvIn, kdIn, betaIn, m)
% Plot the eigenvalue of the R^3 contraction matrix as a function of the
% eigenvalues of the Hessian
% INPUT:
%   eMin := min eigenvalue of the Hessian of the cost function
%   eMax := max eigenvalue of the Hessian of the cost function
%   kvIn := scalar velocity gain
%   kdIn := scalar position gain
%   betaIn := scalar for convergence rate
%   m := the metric tensor as a [2x2] matrix (pos. def)
% OUTPUTS:
%   h := figure handler
h=figure;
Gmat = @(eVal) [2*m(1,2)+betaIn*m(2,2)-2*kvIn*m(2,2), betaIn*m(1,2)-kvIn*m(1,2)-eVal*kdIn*m(2,2)+m(1,1);...
betaIn*m(1,2)-kvIn*m(1,2)-eVal*kdIn*m(2,2)+m(1,1), betaIn*m(1,1)-2*eVal*kdIn*m(1,2)];
hold on
for i = eMin:.001:eMax
    G = Gmat(i);
    eG = sort(eig(G),1,'descend');
    if ( all (eG<=0) )
        plot3(i,eG(1),eG(2),'gx');
    else
        plot3(i,eG(1),eG(2),'rx');
    end
end
xlabel('S'); % eigenvalue of the Hessian
ylabel('eVal 1'); % an eigenvalue of the contracton matrix
zlabel('eVal 2'); % an eigenvalue of the contracton matrix
end

