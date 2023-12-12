%Pack the state of a rotation and its velocity into a vector
%function x=rotDyn_statePack(R,w)
%The vectors x are divided as follows:
%   x(1:9)      Rotation R ([3 x 3] matrix)
%   x(10:12)    Angular velocity w in local coordinates ([3 x 1] vector)
function x=rotDyn_statePack(R,w,varargin)
flagAugmentedSystem = false;
% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'augmentedsystem'
            % The dyanmics is on TSO(3)xSO(3). RReference should be read from x(13:end)
            ivarargin=ivarargin+1;
            RRef = varargin{ivarargin};
            flagAugmentedSystem = true;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

NR=size(R,3);
x=zeros(12,NR);
x(1:9,:)=reshape(R,9,[]);
x(10:12,:)=w;

if flagAugmentedSystem
    % add entries 13:21 for the reference rotation
    x(13:21,:) =reshape(RRef,9,[]);
end