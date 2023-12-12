function [R,w,varargout]=rotDyn_stateUnpack(x,varargin)
flagAugmentedSystem = false;
% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'augmentedsystem'
            % The dyanmics is on TSO(3)xSO(3). RReference should be read from x(13:end)
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
R=reshape(x(1:9,:),3,3,[]);
w=x(10:12,:);

if flagAugmentedSystem
    % Extract the reference rotation from index 13:21
    RRef = reshape(x(13:21,:),3,3,[]);
    varargout{1} = RRef;
end