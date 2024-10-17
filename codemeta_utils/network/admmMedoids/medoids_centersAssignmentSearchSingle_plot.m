function medoids_centersAssignmentSearchSingle_plot(x,mu,kCenter,varargin)
xEval=linspace(0,1,100);
muEval=@(xEval) [mu(:,1:(kCenter-1)) xEval mu(:,(kCenter+1):end)];
f=@(xEval) medoids_centersAssignmentCost(x,muEval(xEval),varargin{:});
funImage(xEval,xEval,f,'method','surf');
