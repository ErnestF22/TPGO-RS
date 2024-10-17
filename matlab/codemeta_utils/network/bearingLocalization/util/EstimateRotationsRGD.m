function  [Rb,cost] = EstimateRotationsRGD(Adj,F,M,maxIter)

flagDebug = false;

if nargin < 4
    maxIter =  2000;
end


sigma_r = 0.01;

% parameters needed: Adj, maxIter, step 
nV = size(Adj,1);

symAdj = (.5*(Adj + Adj') + 10^(-5)) > 1;   

nE = sum(symAdj(:))/2;
  
nF = zeros(numel(F),1);

for iF=1:numel(F)
    nF(iF) = numel(F{iF});
end
dmax = max(sum(symAdj));
kmax =max(nF);
sum_nF = sum(nF);

maxstep = 0.5/((1+kmax)*dmax); 


% random initialization
Rb = zeros(3,3,nV);
for iV= 1:nV
   Rb(:,:,iV) = randrot(3); 
end

% step size definition
step = 0.9*maxstep;



cost = zeros(maxIter,1);

for iIter=1:maxIter
   
    cost(iIter) = .5*trace(myreshape3d2dcolumnwise(Rb)'*M*myreshape3d2dcolumnwise(Rb)) + nE + sum_nF;
    
    my_egrad = 2*myreshape3d2dcolumnwise_inv( M*myreshape3d2dcolumnwise(Rb) );
    
    for iV= 1:nV
        
        rgrad_iV = Rb(:,:,iV)*multiskew(Rb(:,:,iV)'*my_egrad(:,:,iV) + sigma_r*randn(3));
        
        Rb(:,:,iV) = Rb(:,:,iV) * expm( - step*Rb(:,:,iV)'*rgrad_iV   );
        
    end
end

Rb = multitransp(Rb); % don't forget that

if flagDebug
    figure,plot(cost)
end



% some debugging here

mymanifold = rotationsfactory(3, nV);
problem.M = mymanifold;

problem.cost = @(X)  .5*trace(myreshape3d2dcolumnwise(X)'*M*myreshape3d2dcolumnwise(X))+ nE + sum_nF;
problem.egrad = @(X)   myreshape3d2dcolumnwise_inv( M*myreshape3d2dcolumnwise(X) )  ; 
problem.ehess = @(X,U)   myreshape3d2dcolumnwise_inv( M*myreshape3d2dcolumnwise(mymanifold .tangent2ambient(X,U)) )  ; 

 
if flagDebug 
    display('Eigenvalues of the Hessian for gradient descent solution')
    lambdas = hessianspectrum(problem, multitransp(Rb))

end


end

