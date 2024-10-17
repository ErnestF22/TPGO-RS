%function ho_matrix(fid,matrix,varargin)
%Write a matrix as a table
%Arguments:
%   fid     file handle to write to
%   matrix  matrix to write
%Optional arguments:
%   'format',format     format to use for the fprintf
%
%See also ho_tablestart, ho_tablerow, ho_tableend

%%AUTORIGHTS%%

function ho_matrix(fid,matrix,varargin)
format='%.3e';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'format'
            ivarargin=ivarargin+1;
            format=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

ho_tablestart(fid)
for row=1:size(matrix,1)
    ho_tablerow(fid,[],matrix(row,:),'format',format)
end
ho_tableend(fid)