function [X] = solve_manopt(Adj,F,M)

flagDebug = false;

nV = size(M,1)/3;

symAdj = (.5*(Adj + Adj') + 10^(-5)) > 1;   

nE = sum(symAdj(:))/2;
  
nF = zeros(numel(F),1);

for iF=1:numel(F)
    nF(iF) = numel(F{iF});
end

sum_nF = sum(nF);


mymanifold = rotationsfactory(3, nV);
problem.M = mymanifold;

problem.cost = @(X)  .5*trace(myreshape3d2dcolumnwise(X)'*M*myreshape3d2dcolumnwise(X))+ nE + sum_nF;;
problem.egrad = @(X)   myreshape3d2dcolumnwise_inv( M*myreshape3d2dcolumnwise(X) )  ; 
problem.ehess = @(X,U)   myreshape3d2dcolumnwise_inv( M*myreshape3d2dcolumnwise(mymanifold .tangent2ambient(X,U)) )  ; 

if flagDebug
    %checkgradient(problem); pause(0.5);
    %checkhessian(problem); pause(0.5);
end



% Solve the problem
options.tolgradnorm = 10^(-6);
options.verbosity = 0;
X = trustregions(problem,[],options);

if flagDebug
    display('Eigenvalues of the Hessian for manopt solution')
    lambdas = hessianspectrum(problem, X)
end


X = multitransp(X);

end

