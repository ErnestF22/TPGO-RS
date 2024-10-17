%function ho_figure(fid,name,varargin)
%Put an image in an HTML file
%Arguments:
%   fid     file handle to write to
%   name    name of the file with the figure (without extension)
%Optional arguments:
%   'ext',ext       file extension. Default: 'png'
%   'path',path     path to the file as will appear inside the HTML file.
%                   Default: none (directory of the HTML file)
%   'width',width   Specify "width" field in HTML tag. Default: none
%   'currentFigure' Use savefigure to save the current figure in addition
%                   to inserting the HTML
%
%See also savefigure
%

%%AUTORIGHTS%%

function ho_figure(fid,name,varargin)
width=[];
path='';
ext='png';
flagSaveCurrentFigure=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'width'
            ivarargin=ivarargin+1;
            width=varargin{ivarargin};
        case 'path'
            ivarargin=ivarargin+1;
            path=varargin{ivarargin};
        case 'ext'
            ivarargin=ivarargin+1;
            ext=varargin{ivarargin};
        case 'currentfigure'
            ivarargin=ivarargin+1;
            flagSaveCurrentFigure=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagSaveCurrentFigure
    savefigure(name,ext);
end

fprintf(fid,'<p><img ');
if ~isempty(width)
    fprintf(fid, 'width=%dpx ', width);
end
fprintf(fid,'src="%s" /></p>\n',fullfile(path,[name '.' ext]));
