%function ho_tableend(fid)
%End a table in an HTML file
%Arguments:
%   fid     file handle to write to

%%AUTORIGHTS%%

function ho_tableend(fid)
fprintf(fid,'</table>\n');

