%Numerically check the derivative of a function
%function [varargout]=check_der(F,dF,t,varargin)
%Inputs
%   F   handle to function to evaluate
%   dF  handle to expected derivative of f. If the string 'function' is
%       used, the function F should have prototype [f,df]=F(t), where df is
%       the value for dF
function [varargout]=funCheckDer(F,dF,t,varargin)
flagDisplayError=true;
derPlotOpts={};

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(char(lower(varargin{ivarargin})))
        case 'nodisplay'
            flagDisplayError=false;
        case 'derplotopts'
            ivarargin=ivarargin+1;
            derPlotOpts=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

% If it does NOT exist a variable named 't' create it 
if ~exist('t','var')
    % Create a linearly spaced vector t between -1 and 1 of size 100
    t=linspace(-1,1,100);
end

% If it does NOT exist a variable named 'dF' create it 
if ~exist('dF','var') || isempty(dF)
    dF='function';
end

flagDerFromFunction=false;
if ischar(dF)
    switch dF
        case 'function'
            flagDerFromFunction=true;
        otherwise
            error('String for dF not recognized')
    end
end

if ischar(t)
    switch t
        case 'angle'
            t=linspace(-pi,pi,101);
        otherwise
            error('String for t not recognized')
    end
end

d=size(F(0));
if(prod(d)==1)
    Ft=zeros(size(t));
    dFt=zeros(size(t));
    for p=1:length(t)
        if flagDerFromFunction
            [Ft(p),dFt(p)]=F(t(p));
        else
            Ft(p)=F(t(p));
            dFt(p)=dF(t(p));

        end
    end
    appder=funApproxDer(F,t);

    if flagDisplayError
        disp(['Maximum absolute error: ' num2str(max(abs(appder-dFt)))])
        disp(['Maximum relative error: ' num2str(max(abs((appder-dFt)./appder)))])
    end

    flagHold=ishold();
    plot(t,Ft)
    hold on
    plot(t,dFt,'r',t,appder,'x',derPlotOpts{:});
    legend('fun','der','approx der')
    if ~flagHold
        hold off
    end
    grid on
else
    Nt=length(t);
    Ft=zeros(d(1),d(2),Nt);
    dFt=zeros(d(1),d(2),Nt);
    appder=zeros(d(1),d(2),Nt);
    for id=1:prod(d)
        subplot(d(1),d(2),id)
        [iResult,jResult]=ind2sub(d,id);
        ei=reshape([zeros(id-1,1);1;zeros(prod(d)-id,1)],d);
        if flagDerFromFunction
            [Ft(iResult,jResult,:),dFt(iResult,jResult,:),appder(iResult,jResult,:)]=...
                funCheckDer(@(t) eiFt(ei,F,t),'function',t,varargin{:});
        else
            eiF=@(t) multitrace(multiprod(multitransp(ei), F(t)));
            eidF=@(t) multitrace(multiprod(multitransp(ei), dF(t)));
%             eiF=@(t) trace(ei'*F(t));
%             eidF=@(t) trace(ei'*dF(t));
            [Ft(iResult,jResult,:),dFt(iResult,jResult,:),appder(iResult,jResult,:)]=...
                funCheckDer(eiF,eidF,t,varargin{:});
        end

    end
end

if nargout>=1
    varargout{1}=Ft;
    if nargout>=2
        varargout{2}=dFt;
        if nargout>=3
            varargout{3}=appder;
        end
    end
end

function [eiF,eidF]=eiFt(ei,F,t)
[Ft,dFt]=F(t);
eiF=trace(ei'*Ft);
eidF=trace(ei'*dFt);
