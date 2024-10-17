%function ho_tablerow(fid,header,datas,varargin)
%Write a single row in an HTML table
%Arguments:
%   fid     file handle to write to
%   header  header for the row (if empty, no header)
%   datas   array or cell array with the data
%Optional arguments:
%   'format',format     format string to use with fprintf. Default: '%.3f'
%                       if datas is a numeric array, '%s' otherwise

%%AUTORIGHTS%%

function ho_tablerow(fid,header,datas,varargin)
if isnumeric(datas)
    format='%.3f';
else
    format='%s';
end

if isnumeric(datas)
    datas=num2cell(datas);
elseif ischar(datas)
    datas=cellstr(datas);
end

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


fprintf(fid,'<tr>');
if ~isempty(header)
    fprintf(fid,'<th>%s</th>',header);
end

for i=1:length(datas)
    fprintf(fid,['<td align="right">' format '</td>'], datas{i});
end
fprintf(fid,'</tr>\n');
