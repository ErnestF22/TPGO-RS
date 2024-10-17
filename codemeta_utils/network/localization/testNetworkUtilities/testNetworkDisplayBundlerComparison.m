function testNetworkDisplayBundlerComparison(t_node,varargin)
flagFlipZTruth=true;
cameraScale=1;
optsDisplayEst={'Member','gi','Estimated'};
optsDisplayTruth={'DisplayEdges'};
optsDisplayLines={};
methodAbsolutePoses='reference';
flagDisplayCorrespondence=false;
flagDisplayNetwork=false;
displayNetworkOffset=[];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=varargin{ivarargin};
        case 'optsdisplaytruth'
            ivarargin=ivarargin+1;
            optsDisplayTruth=[optsDisplayTruth varargin{ivarargin}];
        case 'optsdisplayest'
            ivarargin=ivarargin+1;
            optsDisplayEst=[optsDisplayEst varargin{ivarargin}];
        case 'optsdisplaylines'
            ivarargin=ivarargin+1;
            optsDisplayLines=[optsDisplayLines varargin{ivarargin}];
        case 'flagflipztruth'
            ivarargin=ivarargin+1;
            flagFlipZTruth=varargin{ivarargin};
        case 'camerascale'
            ivarargin=ivarargin+1;
            cameraScale=varargin{ivarargin};
        case 'displaycorrespondence'
            flagDisplayCorrespondence=true;
            optsDisplayEst=[optsDisplayEst 'flagDisplayEdges' false];
            optsDisplayTruth=[optsDisplayTruth 'flagDisplayEdges' false];
        case 'displaynetworkoffset'
            ivarargin=ivarargin+1;
            displayNetworkOffset=varargin{ivarargin};
            flagDisplayNetwork=true;
        case 'points'
            ivarargin=ivarargin+1;
            X=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

optsDisplayCommon={'OptionsDrawCamera',{'scale',cameraScale,'flagAxes',false},...
    'methodAbsolutePoses',methodAbsolutePoses};
optsDisplayEst=[optsDisplayEst {'OptionsDrawCamera',{'FlipZ'}}];
if flagFlipZTruth
    optsDisplayTruth=[optsDisplayTruth {'OptionsDrawCamera',{'FlipZ'}}];
end    
t_node=testNetworkCompensate(t_node,'methodAbsolutePoses',methodAbsolutePoses);
set(gcf,'Name','Network')
if exist('X','var')
    optsDisplayTruth=[optsDisplayTruth ,'Points',X];
end
testNetworkDisplay(t_node,optsDisplayCommon{:},optsDisplayTruth{:});
hold on
testNetworkDisplay(t_node,optsDisplayCommon{:},optsDisplayEst{:});
if flagDisplayCorrespondence || flagDisplayNetwork
    E=testNetworkGetEdges(t_node);
    switch methodAbsolutePoses
        case 'reference'
            TEst=G2T((cat(3,t_node.gi)));
            TTruth=G2T((cat(3,t_node.gitruth)));
        case 'pose'
            TEst=G2T(invg(cat(3,t_node.gi)));
            TTruth=G2T(invg(cat(3,t_node.gitruth)));
    end
end
if flagDisplayCorrespondence
    plot3([TEst(1,:);TTruth(1,:)],[TEst(2,:);TTruth(2,:)],[TEst(3,:);TTruth(3,:)],'Color',[1 0.7 0],optsDisplayLines{:})
end
if flagDisplayNetwork
    plot3([TTruth(1,E(:,1));TTruth(1,E(:,2))]+displayNetworkOffset(1),...
        [TTruth(2,E(:,1));TTruth(2,E(:,2))]+displayNetworkOffset(2),...
        [TTruth(3,E(:,1));TTruth(3,E(:,2))]+displayNetworkOffset(3),'b:',optsDisplayLines{:})
end    

hold off
axis equal
