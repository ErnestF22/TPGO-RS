function varargout=tagsDatasetError(t_node,varargin)
flagDisplay=false;
color='r';
NTags=6;
tolVarianceTTruth=1e-12;

if nargout==0
    flagDisplay=true;
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'color'
            ivarargin=ivarargin+1;
            color=varargin{ivarargin};
        case 'ntags'
            ivarargin=ivarargin+1;
            NTags=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

TOpt=subMatrix(G2T(cat(3,t_node.gi)),1:3,1:NTags);
TTruth=subMatrix(G2T(cat(3,t_node.gitruth)),1:3,1:NTags);
ROpt=G2R(cat(3,t_node.gi));
ROpt=ROpt(:,:,1:NTags);
RTruth=G2R(cat(3,t_node.gitruth));
RTruth=RTruth(:,:,1:NTags);

flagDegenerateCase=trace(cov(TTruth))<tolVarianceTTruth;

if flagDegenerateCase
    disp('Degenerate configutation')
    TOptTransf=TOpt-(mean(TOpt,2)+mean(TTruth,2))*ones(1,NTags);
    
    dSqrt=sqrt(sum((TOptTransf(:)-TTruth(:)).^2));
    a=mean(rotationProcrustesAlignError(ROpt,RTruth));
else
    [d,TOptTransf,transf]=procrustes(TTruth',TOpt','scaling',false,'reflection',false);
    TOptTransf=TOptTransf';
    dSqrt=sqrt(d);
    a=mean(rot_dist(RTruth,multRotArray(transf.T',ROpt),'vector'))*180/pi;
end

if flagDisplay
    disp('Sqrt procrustes error [m]')
    disp(dSqrt)

    plotPoints(TTruth)
    hold on
    plotPoints(TOptTransf,{'Color',color})
    hold off
    axis equal
    
    disp('Mean rotation error [deg]')
    disp(a)
end

if nargout>0
    varargout{1}=dSqrt;
    if nargout>1
        varargout{2}=a;
    end
end
