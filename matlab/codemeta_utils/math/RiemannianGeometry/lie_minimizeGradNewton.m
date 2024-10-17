function [y,output]=lie_minimizeGradNewton(lf,f,y0,varargin)
flagGetErrors=false;        %record stats in output variable
flagDisplayIt=false;
flagProgressBar=false;

flagPhase12=true;

methodShow='none';
stopCriterion='normgrad';
stopThreshold=1e-12;

maxIt=100;
minIt=1;
phase2NPointMax=2;
phase1ResetCounterMax=5;
phase3IncreaseCounterMax=10;

rho=0.8;
epsilon=0.1;
tolMinLineSearch=1e-12;

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
        case 'geterrors'
            flagGetErrors=true;
        case 'showcost'
            methodShow='lin';
        case 'showlogcost'
            methodShow='log';
         case 'stopcriterion'
            ivarargin=ivarargin+1;
            stopCriterion=lower(varargin{ivarargin});
         case 'stopthreshold'
            ivarargin=ivarargin+1;
            stopThreshold=lower(varargin{ivarargin});
        case 'rho'
            ivarargin=ivarargin+1;
            rho=lower(varargin{ivarargin});
        case 'epsilon'
            ivarargin=ivarargin+1;
            epsilon=lower(varargin{ivarargin});
        case {'disablephase12','gradientonly'}
            flagPhase12=false;
        case 'progressbar'
            flagProgressBar=true;
        case 'flagprogressbar'
            ivarargin=ivarargin+1;
            flagProgressBar=lower(varargin{ivarargin});
        case 'phase3increasecountermax'
            ivarargin=ivarargin+1;
            phase3IncreaseCounterMax=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagDisplayIt
    flagProgressBar=false;
end

if ~strcmpi(methodShow,'none')
    flagGetErrors=true;
end

if flagGetErrors
    allc=NaN(maxIt,1);
    allt=NaN(maxIt,1);
end

y=y0;

t0=cputime;
if flagPhase12
    [c,g,H]=f(y);
    delta=computeDelta(H);
else
    [c,g]=f(y);
    delta=NaN;
end
szGradient=size(g);
phase1ResetCounter=0;
phase3IncreaseCounter=0;
flagRunning=true;
if flagProgressBar
    w=getTextWaitBar(maxIt);
end
for it=1:maxIt
    if flagDisplayIt
        deltaAction='';
        ng=norm(g);
    end
    if flagDisplayIt || flagGetErrors
        tCurrent=cputime-t0;
    end
    if flagGetErrors
        allt(it)=tCurrent;
        allc(it)=c;
    end
    
    flagPhase1Success=false;
    flagPhase2Success=false;
    if flagPhase12
        %compute Levemberg-Marquardt direction
        D=H+delta*eye(size(H));
        d=reshape(-D\g(:),szGradient);

        %check new cost
        yNew=lf.exp(y,lf.hat(y,d));
        [cNew,gNew,HNew]=f(yNew);
        if cNew<c
            %LM was a success
            phase1ResetCounter=0;
            flagPhase1Success=true;
            delta=delta*rho;
            if flagDisplayIt
                deltaAction=[deltaAction 'v'];
            end
        else
            phase1ResetCounter=phase1ResetCounter+1;
            if phase1ResetCounter==phase1ResetCounterMax
                %increased delta too many times, just reset it
                phase1ResetCounter=0;
                delta=computeDelta(H);
                if flagDisplayIt
                    deltaAction=[deltaAction 'R'];
                end
            else
                delta=delta/rho;
                if flagDisplayIt
                    deltaAction=[deltaAction '^'];
                end
            end
        end
        if ~flagPhase1Success
            %phase 1 failed, go to phase 2
            %linear interpolation between LM and gradient directions
            dLM=d;
            flagPhase2Success=false;
            for phase2NPoint=1:phase2NPointMax
                if flagDisplayIt
                    deltaAction=[deltaAction '.'];
                end
                t=phase2NPoint/(phase2NPointMax+1);
                d=(1-t)*dLM-t*epsilon*g;
                yNew=lf.exp(y,lf.hat(y,d));
                [cNew,gNew,HNew]=f(yNew);
                if cNew<c
                    flagPhase2Success=true;
                    break;
                end
            end
        end
    end
    if ~flagPhase1Success && ~flagPhase2Success
        %phase 2 failed, go to phase 3 (gradient descent)
        flagPhase3Success=false;
        flagPhase3FirstTrial=true;
        ng=norm(g);
        if phase3IncreaseCounter==phase3IncreaseCounterMax
            if flagDisplayIt
                deltaAction=[deltaAction '~'];
            end
            phase3IncreaseCounter=0;
            epsilon=epsilon/rho;
        end            
        while epsilon*ng>tolMinLineSearch || flagPhase3FirstTrial
            if flagDisplayIt
                deltaAction=[deltaAction '_'];
            end
            d=-epsilon*g;
            yNew=lf.exp(y,lf.hat(y,d));
            if flagPhase12
                [cNew,gNew,HNew]=f(yNew);
            else
                [cNew,gNew]=f(yNew);
            end                
            if cNew<c
                flagPhase3Success=true;
                if flagPhase3FirstTrial
                    phase3IncreaseCounter=phase3IncreaseCounter+1;
                else
                    phase3IncreaseCounter=0;
                end
                break;
            end
            flagPhase3FirstTrial=false;
            epsilon=rho*epsilon;
        end
        if ~flagPhase3Success
            %step size too small, give up
            flagRunning=false;
        end
    end
    if flagDisplayIt
        fprintf('it=%4d t=%5.1fs c=%.4e ng=%.4e delta=%.4e epsilon=%.4e (%s)\n',it,tCurrent,c,ng,delta,epsilon,deltaAction);
    end
    
    if flagRunning && it>minIt && ~strcmpi(stopCriterion,'none')
        switch stopCriterion
            case 'fun'
                if c-cNew<stopThreshold
                    flagRunning=false;
                end
            case 'normgrad'
                if norm(g)<stopThreshold
                    flagRunning=false;
                end
            otherwise
                error('Stopping criterion not recognized')
        end
    end
                
    y=yNew;
    c=cNew;
    g=gNew;
    if flagPhase12
        H=HNew;
    end
    
    if ~flagRunning
        break
    end
    
    if flagProgressBar
        w(it)
    end
end
if flagProgressBar
    fprintf('\n');
end

if flagGetErrors
    allc=allc(1:it);
    allt=allt(1:it);
    switch methodShow
        case 'log'
            subplot(1,2,1)
            semilogy(allc)
            subplot(1,2,2)
            semilogy(allt,allc)
        case 'lin'
            subplot(1,2,1)
            plot(allc)
            subplot(1,2,2)
            plot(allt,allc)
    end
end

%If requested, prepare output variable with additional info
if nargout>1
    output.cFinal=c;
    output.gFinal=g;
    if flagPhase12
        output.HFinal=H;
    end
    output.itReached=it;
    output.flagMaxItReached=(it==maxIt);
    if flagGetErrors
        output.cost=allc;
        output.t=allt;
    end
end

function delta=computeDelta(H)
delta=-min(min(diag(H)-sum(abs(H-diag(diag(H))),2)),0);
