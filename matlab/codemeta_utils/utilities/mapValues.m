%function vMapped=mapValues(v,m,varargin)
%Maps values in v from m(:,1) to m(:,2)
%If the 
%If m has only one column, the second is filled to be (1:size(m,1))'
%If m is omitted or empty, the first column is filled with unique(v), and
%the second as the above.
%Optional inputs
%   'actionMapped',s    Determines what to do for values in v are not found in m(:,1)
%       'none'          Values are left unchanged (default)
%       'remove'        Values are removed from vMapped

function vMapped=mapValues(v,m1,varargin)
actionUnmapped='none';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'actionunmapped'
            ivarargin=ivarargin+1;
            actionUnmapped=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~exist('m1','var') || isempty(m1)
    m1=shiftdim(unique(v));
end

if size(m1,2)==1
    nm=size(m1);
    m1(:,2)=(1:nm)';
end

vMapped=v;

switch actionUnmapped
    case 'none'
        %nothing to do
    case 'remove'
        vMapped(~ismember(vMapped,m1(:,1)))=[];
    otherwise
        error('actionUnmapped not recognized')
end

for im=1:size(m1,1)
    vMapped(v==m1(im,1))=m1(im,2);
end

