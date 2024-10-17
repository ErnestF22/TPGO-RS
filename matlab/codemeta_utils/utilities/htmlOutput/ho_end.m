%function ho_end(fid)
%Write the footer of an HTML file
%Arguments:
%   fid     file handle to write to

%%AUTORIGHTS%%

function ho_end(fid)
fprintf(fid,'</body>\n</html>\n');
