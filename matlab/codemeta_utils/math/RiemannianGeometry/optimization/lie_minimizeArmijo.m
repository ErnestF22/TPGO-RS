function [x,errors]=lie_minimizeArmijo(lf,f,gradf,x0,varargin)

    DEBUG=false;     %activate debug checks (e.g. orthogonality)
    
    itersMax=2000;  %max number of iterations
    tolGrad=1e-14;  %minimum norm of the gradient before stopping
    tolFun=1e-14;   %minimum function decrease before stopping
    
    armijoOpts={};  %options for Armijo search
    retractionsOpts={}; %options for computing retractions
    
    methodDirection='conjgrad';  %method to compute search direction
    methodLineSearch='armijo';   %method to conduct line search
    
    %init
    x=x0;
    feval=f(x);
    dim=lf.dim(x0);   %dimension of the manifold
    
    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'itersmax'
                ivarargin=ivarargin+1;
                itersMax=varargin{ivarargin};
            case 'armijoopts'
                ivarargin=ivarargin+1;
                armijoOpts=varargin{ivarargin};
            case 'retractionsopts'
                ivarargin=ivarargin+1;
                retractionsOpts=varargin{ivarargin};
            case 'methoddirection'
                ivarargin=ivarargin+1;
                methodDirection=varargin{ivarargin};
            case 'methodlinesearch'
                ivarargin=ivarargin+1;
                methodLineSearch=varargin{ivarargin};
            case 'debug'
                DEBUG=true;
            otherwise
                disp(varargin{ivarargin})
                error('Argument not valid!')
        end
        ivarargin=ivarargin+1;
    end
    
    %prepare vectors to store errors
    stepsize=Inf(1,itersMax);
    allfeval=Inf(1,itersMax);
    allm=Inf(1,itersMax);
    alltimes=Inf(1,itersMax);
    allNormGrad=Inf(1,itersMax);

    tstart=cputime;
    %iterations
    for it=1:itersMax
        allfeval(it)=feval;
        alltimes(it)=cputime-tstart;
        
        %gradient
        gf=gradf(x);
        normgf=sqrt(lf.metric(x,gf,gf));
        allNormGrad(it)=normgf;
        
        if DEBUG
            gfProj=lf.tangentProj(x,gf);
            avgErr=sum(sum(abs(gf-gfProj)))/numel(gf);
            if avgErr>1e-11
                warning(['Gradient does not appear to be a tangent vector (avgErr ' num2str(avgErr) ')'])
                keyboard
            end
        end

        
        if normgf<tolGrad
            errors.exitInfo=['Minimum gradient norm reached: ' num2str(normgf)];
            break;
        end
            
        
        %descent direction
        switch methodDirection
            case 'gradient'
                eta=-gf;
            case 'conjgrad'         %conjugate gradient with approximated parallel transport
                if mod(it,dim)==1
                    eta=-gf;
                else
                    %gfprev=lf.tangentProj(x,gfprev);
                    %eta=lf.tangentProj(x,eta);
                    
                    gamma=lf.metric(x,gf-0,gf)/lf.metric(xprev,gfprev,gfprev);  %here use parallel transport of gk=0
                    neweta=-gf+lf.tangentProj(x,gamma*eta);
                    if lf.metric(x,gf,neweta)<0     %check that we got a descent direction
                        eta=neweta;
                    else
                        eta=-gf;
                        disp(['Reset direction at iteration ' num2str(it)])
                    end
                end
                gfprev=gf;
        end
        
        fprev=feval;
        xprev=x;
        switch methodLineSearch
            case 'armijo'
                %Armijo inexact line search
                [x,feval,t,m,failFlag,infoArmijoSearch]=lie_armijoSearch(lf,f,x,feval,gf,eta,retractionsOpts,armijoOpts{:});
                allm(it)=m;
                stepsize(it)=t*normgf;

                if DEBUG
                    dbCheckOrthogonality(x);
                end

                %Checks
                if failFlag
                    errors.exitInfo=['Couldn''t find Armijo point during line search'];
                    errors.infoArmijoSearch=infoArmijoSearch;
                    break
                end
            case 'fixedStepsize'
                x=lf.retractions(x0,eta,retractionsOpts{:});
                feval=f(x);
            case 'exact'
                %find point with greater cost than current
                etaNext=eta;
                fnext=f(lf.retractions(x0,etaNext,retractionsOpts{:}));
                while fnext<feval+1e-13
                    etaNext=2*etaNext;
                    fnext=f(lf.retractions(x0,etaNext,retractionsOpts{:}));
                end
                %exact minimization along the geodesic
                fretrac=@(t) f(lf.retractions(x0,t*etaNext,retractionsOpts{:}));
                [t,feval,exitFlag,infoSearch] = fminbnd(fretrac,0,1);
                x=lf.retractions(x0,t*etaNext,retractionsOpts{:});
        end
                
        if fprev<feval
            warning(['Cost went up at iteration ' num2str(it)])
            disp(['fprev-feval=' num2str(feval-fprev)])
            check_RiemGrad(lf,f,gradf,xprev);
            disp(['-<gf,eta>=' num2str(-lf.metric(xprev,gf,eta))])
            
            fretrac=@(t) f(lf.retractions(xprev,t*eta,retractionsOpts{:}));
            %fflat=@(t) f(xprev+t*eta);
            
            figure(1000)
            plotfun(fretrac,linspace(-0.1,0.1,100),'.-')
            keyboard
        end
        
        if fprev-feval<tolFun
            errors.exitInfo=['Minimum function decrease reached: ' num2str(fprev-feval)];
            break
        end

    end

    if it==itersMax
        errors.exitInfo=['Maximum number of iterations ' num2str(itersMax) ' reached'];
    end
    
    %prepare errors structure
    errors.stepsize=stepsize(1:it);
    errors.armijoExponent=allm(1:it);
    errors.feval=allfeval(1:it);
    errors.t=alltimes(1:it);
    errors.itNum=it;
    
    %disp(errors.exitInfo)
end



