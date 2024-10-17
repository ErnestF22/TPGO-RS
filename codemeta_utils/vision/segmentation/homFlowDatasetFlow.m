%Generate random infinitesimal flow of given structure
%function [x,dx,v,w]=homFlowDatasetFlow(X,G,varargin)
%Note: v,w are in the body-fixed frame
function [x,dx,v,w]=homFlowDatasetFlow(X,G,varargin)
w=rand(3,1);
v=rand(3,1);
flagNormalizeV=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'v'
            ivarargin=ivarargin+1;
            v=varargin{ivarargin};
        case 'w'
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
if flagNormalizeV
    v=cnormalize(v);
end

[x,JRT]=projectFromG(G,X,'references');


dx=squeeze(multiprod(JRT,[w;v]));
