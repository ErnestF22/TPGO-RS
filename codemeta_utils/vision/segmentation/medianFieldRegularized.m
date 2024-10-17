%Regularized median filtering of arrays
%function uMedian=medianFieldRegularized(uReg,idxList,lambda,varargin)
%Updates each entry with the regularized median of the neighboring entires
%and the regularizer value
%Inputs
%   uReg        value of u for the regularizer
%   idxList     list of linear indeces of the neighbors
%   kList       number of neighbors for each location
%   lambda      cost parameters
%       lambda(1)   for regularizer cost
%       lambda(2)   for neighbors cost
%Optional Inputs
%   maxIt       number of iterations
%   collect     uMedian is a [maxIt x numel(uReg)] array with the value of
%               the estimate at every iteration 
%   rho         overrelaxation parameters (i.e., the stepsize, default: 0.6)
function uMedian=medianFieldRegularized(uReg,idxList,kList,lambda,varargin)
maxIt=10;
flagCollectIterations=false;
rho=0.5;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'collect'
            flagCollectIterations=true;
        case 'rho'
            ivarargin=ivarargin+1;
            rho=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%p is a linear index to the current slice of the output
%If flagCollectIteration is false, p will always be set to zero
%If flagCollectIteration is true, p will be increased through all the
%slices of the output to collect all the iterations
p=0;
Nu=numel(uReg);
if ~flagCollectIterations
    uMedian=uReg;
else
    uMedian=zeros(Nu,maxIt+1);
    uMedian(1:Nu)=uReg;
end
uMedianNext=zeros(size(uReg));
for it=1:maxIt
    for iu=1:Nu
        uMedianNext(iu)=medianRegularized(uMedian(p*Nu+idxList(iu,1:kList(iu))),lambda(2),uReg(iu),lambda(1));
    end
    
    uMedianPrev=uMedian(p*Nu+1:(p+1)*Nu);
    if flagCollectIterations
        p=p+1;
    end
    uMedian(p*Nu+1:(p+1)*Nu)=(1-rho)*uMedianPrev+rho*uMedianNext(:)';
end
