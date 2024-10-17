%function ho_start(fid,title)
%Write the header of an HTML file
%Arguments:
%   fid     file handle to write to
%   title   string of text for webpage title (appears in browser window/tab
%           title. Default: empty

%%AUTORIGHTS%%

function ho_start(fid,title)
fprintf(fid,'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n<html>\n<head>\n');
if(exist('title','var'))
    fprintf(fid,'<title>%s</title>\n',title);
end
fprintf(fid,'</head>\n<body>\n');
