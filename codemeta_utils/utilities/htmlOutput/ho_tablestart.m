%function ho_tablestart(fid,titles,varargin)
%Start a table in an HTML file
%Arguments:
%   fid     file handle to write to
%   titles  if not empty, insert a row with titles. Default: empty
%Optional arguments:
%   'summary',summary   content for 'summary' option of HTML tag

%%AUTORIGHTS%%

function ho_tablestart(fid,titles,varargin)
if ~exist('titles','var')
    titles=[];
end
summary=[];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'summary'
            ivarargin=ivarargin+1;
            summary=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ischar(titles)
    titles=cellstr(titles);
end

fprintf(fid,'<table border="1" cellpadding="5" ');
if ~isempty(summary)
    fprintf(fid,'summary="%s"',summary);
end

fprintf(fid,'>\n');

if ~isempty(titles)
    fprintf(fid,'<tr>');
    for i=1:length(titles)
        fprintf(fid,'<th>%s</th>',titles{i});
    end
    fprintf(fid,'</tr>\n');
end
