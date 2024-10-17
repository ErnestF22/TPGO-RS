%function [yMean,output]=lie_minimize(lf,fun,der,y0,varargin)
%Implements gradient descent on manifolds with Armijo's line search.
%Arguments:
%   lf      struct with function to implement exp map and metric
%   fun     function to be minimized. Must have prototype c=fun(y). Note: y
%           can be a 3-D array for multivariate minimization
%   der     gradient of f. Must have prototype h=der(y), where h is a
%           tangent vector at y. Note: if y is a 3-D array (multivariate
%           minimization), then h must be a 3-D array too with opportune
%           dimensions
%   y0      point in the manifold where to start the minimization
%
%Optional arguments:
%   'maxIt',NIt         Maximum number of iterations
%   'minIt',NIt         Minimum number of iterations
%   'showCost'          Show a linear plot with the cost at every iteration
%                       before terminating
%   'showLogCost'       Same as showCost, but use semilogy for plotting
%   'stepSize',epsilon  Set initial step size to use for line search to
%                       epsilon
%   'lineSearch',flag   If true, perform a line search to find a suitable
%                       step size (this is the default). If false, make a
%                       step of size epsilon.
%   'tolGrad',tol       Tolerance on the norm of the gradient to stop.
%   'tolCost',tol       Tolerance on the cost decrease to stop.
%   'stopThreshold',thr Threshold on cost to stop
%
%   'stopCriterion',str Determines the criterion for stopping the algorithm
%       'standard'          Stop when decrease in cost is less than tolCost
%                           or when norm of gradient is less than tolGrad
%       'nostop'            Run all maxIt iterations
%       'threshold'         Stop when cost falls below stopThreshold     

%%AUTORIGHTS%%

function [yMean,output]=lie_minimize(lf,fun,der,yMean,varargin)
    maxIt=100;      %maximum number of iterations
    minIt=1;        %minimum number of iterations before stopping
    methodShow='none';
    flagLineSearch=true;
    flagDisplayIt=false;        %display stats at every iteration
    flagDisplayProgressBar=false;   %display a progress bar
    flagGetErrors=false;        %record in output variable the cost at every iteration
    stopCriterion='standard';
    stopThreshold=-Inf;
    
    tolGrad=1e-8;
    tolCost=1e-14;

    epsilon=1;

    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'maxit'
                ivarargin=ivarargin+1;
                maxIt=varargin{ivarargin};
            case 'minit'
                ivarargin=ivarargin+1;
                minIt=varargin{ivarargin};
            case 'displayit'
                flagDisplayIt=true;
            case 'displaybar'
                flagDisplayProgress='lin';
            case 'showcost'
                methodShow='lin';
            case 'showlogcost'
                methodShow='log';
            case 'stepsize'
                ivarargin=ivarargin+1;
                epsilon=varargin{ivarargin};
            case 'linesearch'
                ivarargin=ivarargin+1;
                flagLineSearch=varargin{ivarargin};
            case 'geterrors'
                flagGetErrors=true;
            case 'tolgrad'
                ivarargin=ivarargin+1;
                tolGrad=varargin{ivarargin};
             case 'tolcost'
                ivarargin=ivarargin+1;
                tolCost=varargin{ivarargin};
             case 'stopcriterion'
                ivarargin=ivarargin+1;
                stopCriterion=lower(varargin{ivarargin});
             case 'stopthreshold'
                ivarargin=ivarargin+1;
                stopThreshold=lower(varargin{ivarargin});
            otherwise
                disp(varargin{ivarargin})
                error('Argument not valid!')
        end
        ivarargin=ivarargin+1;
    end

    fPrev=Inf;
    if ~strcmpi(methodShow,'none') || flagGetErrors
        allf(1)=fun(yMean);
    end

    if flagDisplayProgressBar
        w=getTextWaitBar(maxIt);
        w(0)
    end
    for it=1:maxIt
        if flagDisplayIt
            fprintf('it=%d',it)
        end
        
        yPrev=yMean;
        hy=-der(yMean);
        
        if flagDisplayIt
            fprintf(' gradNorm=%e', sqrt(sum(lf.metric(yMean,hy,hy))))
        end

        if(size(yMean,3)==1)
            geodCost=@(t) fun(lf.exp(yMean,t*hy));
        else
            geodCost=@(t) evalFunGeod(lf,fun,yMean,t*hy);
        end

        if it>minIt && strcmp(stopCriterion,'standard') && sqrt(sum(lf.metric(yMean,hy,hy).^2))<tolGrad
            break
        end
        
        if flagLineSearch
            topt=fminbnd(geodCost,-epsilon/10,epsilon);  %start from slightly less than zero to allow convergence to stepsize zero
            if(topt<-1e-6)
                warning(['Going backward on the descent direction at iteration ' num2str(it)])
            end
        else
            topt=epsilon;
        end
        if flagDisplayIt
            fprintf(' t=%e', topt);
        end
        
        %disp(topt)
        
        if(size(yMean,3)==1)
            yMean=lf.exp(yMean,topt*hy);
        else
            yMean=multiexp(lf,yMean,topt*hy);
        end

        fCurrent=fun(yMean);
        if flagDisplayIt
            fprintf(' f=%e\n', fCurrent);
        end
        if(fCurrent-fPrev>1e-7)
            warning(['Cost went up at iteration ' num2str(it)...
                ' with gradient magnitude ' num2str(sqrt(lf.metric(yPrev,hy,hy))) '! Stopping here'])
            break;
        end
        
%         if(size(yMean,3)==1)
%             d=lf.dist(yPrev,yMean);
%         else
%             d=multidist(lf,yPrev,yMean);
%         end
        
        if it>minIt && ((strcmp(stopCriterion,'standard') && (fPrev-fCurrent)<tolCost) ...
            || (strcmp(stopCriterion, 'threshold') && fCurrent<stopThreshold))
            %disp(['Terminated after ' num2str(it) ' iteration with cost ' num2str(fCurrent,'%.12e')])
            break;
        end
        
        fPrev=fCurrent;
        if ~strcmpi(methodShow,'none') || flagGetErrors
            allf(it+1)=fun(yMean);
        end
        if flagDisplayProgressBar
            w(it)
        end
    end

    %Check if stopping criterion was not met
    if ~strcmp(stopCriterion,'nostop') && it==maxIt
        disp(['  Minimization did not converge after ' num2str(maxIt) ' iterations'])
    end
    
    %If requested, show cost
    if ~strcmpi(methodShow,'none')
        switch methodShow
            case 'lin'
                plot(allf)
            case 'log'
                semilogy(allf)
        end
    end
    
    %If requested, prepare output variable with additional info
    if nargout>1
        output.fFinal=fCurrent;
        output.hFinal=hy;
        output.stepSizeFinal=topt;
        output.itReached=it;
        output.maxItReached=(it==maxIt);
        if flagGetErrors
            output.cost=allf;
        end
    end

end

%evaluate a function along a geodesic
function f=evalFunGeod(lf,fun,y0,h)
    y=multiexp(lf,y0,h);
    f=fun(y);
end

function y=multiexp(lf,y0,h)
    y=zeros(size(y0));
    for n=1:size(y0,3)
        y(:,:,n)=lf.exp(y0(:,:,n),h(:,:,n));
    end
end        

function d=multidist(lf,y1,y2)
    d=0;
    for n=1:size(y1,3)
        d=d+lf.dist(y1(:,:,n),y2(:,:,n));
    end
end        
   