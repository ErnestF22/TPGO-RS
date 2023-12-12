%function testNetworkDisplayErrors(t_node,type,varargin)
%Display empirical cumulative distributions of errors from results in t_node
%Arguments
%   type    Type of error to plot. It can be one or any combination of:
%       'r'     rotation
%       't'     translation
%       'n'     network
%Optional arguments
%   'boxplots'                  Overimpose box plots to cumulative distributions
%   'optscomputeerrors',opts    Cell array of arguments to be passed to
%                               testNetworkComputeErrors
function output=testNetworkDisplayErrors(t_node,type,varargin)
if ~exist('type','var')
    type='RT';
end
numberFiguresFilled=0;

optsComputeError={};
optsNetworkCompensate={};
flagAddBoxPlots=false;
flagTextSummary=true;
flagTextLatex=false;
methodTranslErr='angle';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'methodtranslationerror'
            ivarargin=ivarargin+1;
            methodTranslErr=lower(varargin{ivarargin});
        case 'boxplots'
            flagAddBoxPlots=true;
        case 'latex'
            flagTextSummary=true;
            flagTextLatex=true;
        case 'optscomputeerrors'
            ivarargin=ivarargin+1;
            optsComputeError=[optsComputeError varargin{ivarargin}];        %#ok<AGROW>
            if ~iscell(optsComputeError)
                optsComputeError={optsComputeError};
            end
        case {'references','poses'}
            optsComputeError=[optsComputeError varargin{ivarargin}];
        case 'methodabsoluteposes'
            optsComputeError=[optsComputeError varargin{ivarargin:ivarargin+1}];
            ivarargin=ivarargin+1;
        case 'optsnetworkcompensate'
            ivarargin=ivarargin+1;
            optsNetworkCompensate=[optsNetworkCompensate varargin{ivarargin}];        %#ok<AGROW>
            if ~iscell(optsNetworkCompensate)
                optsNetworkCompensate={optsNetworkCompensate};
            end
       otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


%recursive call if both rotations and translations errors are requested
%TODO: this is not efficient, computes the errors twice. Maybe move to a
%iterative solution
type=unique(lower(type));
if any(type=='r') && any(type=='t')
    if numberFiguresFilled>0
        figure
    end
    subplot(2,1,1)
    title('Rotation errors')
    outputR=testNetworkDisplayErrors(t_node,'r',varargin{:});

    subplot(2,1,2)
    title('Translation errors')
    outputT=testNetworkDisplayErrors(t_node,'t',varargin{:});

    if nargout>0
        output.R=outputR;
        output.T=outputT;
    end
    
    legend('Localization','Measurements')
    numberFiguresFilled=numberFiguresFilled+1;
elseif any(type=='r') || any(type=='t')
    if numberFiguresFilled>0
        figure
    end

    %compute errors on edges
    [rotErr,translErr,scaleRatios,translErrNorm]=testNetworkComputeErrors(t_node,optsComputeError{:});
    [rotErrMs,translErrMs,scaleRatiosMs,translErrNormMs]=testNetworkComputeErrors(t_node,'Measurements',optsComputeError{:});

    switch lower(type)
        case {'r','rot','rotation'}
            err=rotErr;
            errMs=rotErrMs;
            textSummaryTitle='Rotation';
        case {'t','transl','translation'}
            switch methodTranslErr
                case 'angle'
                    err=translErr;
                    errMs=translErrMs;
                case 'norm'
                    err=translErrNorm;
                    errMs=translErrNormMs;
            end
            textSummaryTitle='Translation';
        otherwise
            error(['Type of error '  type ' not valid!'])
    end

    %convert into degrees
    err=err*180/pi;
    errMs=errMs*180/pi;

    %plot
    cumDistPerc(err,'b')
    hold on
    cumDistPerc(errMs,'r')
    hold off
    h1=gca;

    if flagAddBoxPlots
        pos=get(gca,'Position');
        h2=axes();
        boxplot(gca,[errMs err],'orientation','horizontal','colors','rb')
        set(h2,'Visible','off')
        pos(4)=0.8*pos(4);
        set(h2,'Position',pos)
        set(h2,'YTickLabel',{' '})
        hold off
    end

    set(gcf,'CurrentAxes',h1);
    set(gca,'YTick',0:20:100,'YGrid','on','YMinorGrid','on','YMinorTick','on','XMinorTick','on','YLim',[0 100])
    xlabel('Angular error [deg]')
    ylabel('Percentage of edges')

    if flagTextSummary
        if ~flagTextLatex
            fprintf('\n--- %s errors ---\n',textSummaryTitle);
            fprintf('Localization: mean=%.4e, std=%.4e, median=%.4e, min=%.4e, max=%.4e\n', mean(err), std(err), median(err), min(err), max(err));
            fprintf('Measurements: mean=%.4e, std=%.4e, median=%.4e, min=%.4e, max=%.4e\n', mean(errMs), std(errMs), median(errMs), min(err), max(err));
        else
            fprintf('\\multirow{2}{*}{%s} & Measurements & %.3f & %.3f & %.3f \\\\\n ',  textSummaryTitle, mean(errMs), std(errMs), median(errMs));
            fprintf(' & Estimated & %.3f & %.3f & %.3f \\\\\n ', mean(err), std(err), median(err));
        end
    end
    if nargout>0
        output.err=err;
        output.meanErr=mean(err);
        output.stdErr=std(err);
        output.medianErr=median(err);
        output.errMs=errMs;
        output.meanErrMs=mean(errMs);
        output.stdErrMs=std(errMs);
        output.medianErrMs=median(errMs);
    end    
    numberFiguresFilled=numberFiguresFilled+1;
end

if any(type=='s')
    if ~exist('scaleRatios','var')
        [~,~,scaleRatios]=testNetworkComputeErrors(t_node,optsComputeError{:});
        [~,~,scaleRatiosMs]=testNetworkComputeErrors(t_node,'Measurements',optsComputeError{:});
    end
    fprintf('\n--- Scale Errors ---\n')
    fprintf('Localization: geometric median=%.4e\n',geostd(scaleRatios));
    fprintf('Measurements: geometric median=%.4e\n',geostd(scaleRatiosMs));
end    
        
if any(type=='n')
    %if we already plotted errors, open a new figure
    if numberFiguresFilled>0
        figure
    end
    
    %extract interpretation of poses from optsComputeError
    optsNetworkDisplay=optionsExtract(optsComputeError,{'Poses',0,'References',0,'methodabsoluteposes',1});
    flagHasTruth=isfield(t_node,'gitruth');
    if flagHasTruth
        t_node=testNetworkCompensate(t_node,optsNetworkDisplay{:},optsNetworkCompensate{:});%'Poses');
        testNetworkDisplay(t_node,'DisplayEdges',optsNetworkDisplay{:});
        hold on
    end
    testNetworkDisplay(t_node,'Member','gi','Estimated',optsNetworkDisplay{:});
    if flagHasTruth
        hold off
    end
    numberFiguresFilled=numberFiguresFilled+1;
end    

