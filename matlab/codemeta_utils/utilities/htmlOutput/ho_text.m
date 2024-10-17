%function ho_text(fid,text)
%Print a paragraph in an HTML file
%Arguments:
%   fid     file handle to write to
%   text    string of text to insert

%%AUTORIGHTS%%

function ho_text(fid,text)
fprintf(fid,'<p>%s</p>\n',text);
