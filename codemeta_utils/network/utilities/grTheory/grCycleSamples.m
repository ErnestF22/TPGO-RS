function [C,output]=grCycleSamples(E,varargin)
maxReps=100;
NSamples=maxReps;
flagCollect=nargout>1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'collect'
            flagCollect='reference';
        case 'maxreps'
            ivarargin=ivarargin+1;
            maxReps=varargin{ivarargin};
        case 'nsamples'
            ivarargin=ivarargin+1;
            NSamples=varargin{ivarargin};
            if maxReps<NSamples
                maxReps=NSamples;
            end
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

C=grCycleBasis(E);
if flagCollect
    NCycles=zeros(1,NReps+1);
    NCycles(1)=size(C,2);
end
w=sum(C,2);
for iRep=1:maxReps
    CAdd=grCycleBasis([E(:,1:2) w]);
    C=unique([C CAdd]','rows')';
    w=w+sum(CAdd,2);
    if flagCollect
        NCycles(iRep+1)=size(C,2);
    end
    if size(C,2)>NSamples
        C=C(:,1:NSamples);
        if flagCollect
            NCycles=NCycles(1:iRep+1);
        end
        break
    end
end

if flagCollect
    output.NCycles=NCycles;
end
