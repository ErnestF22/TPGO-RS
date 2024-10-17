%function a=getelement(a,stridx)
%Returns eval(['a=a' stridx]);
function a=getelement(a,stridx)
eval(['a=a' stridx ';']);
