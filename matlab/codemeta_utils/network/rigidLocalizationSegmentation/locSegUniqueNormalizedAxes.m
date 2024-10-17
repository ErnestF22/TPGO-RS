%Keep the vectors in v which are unique up to scale and normalize them
%function vUniqueNorm=locSegUniqueNormalizedAxes(v)
function [vUniqueNorm,flagsVUniqueNorm]=locSegUniqueNormalizedAxes(v)
tolCosine=1e-6;

NV=size(v,2);
[c,vNorm]=computeAllCosines(v);
boolDistance=1-abs(c)<tolCosine;
[flagsVUniqueNorm,iv]=unique(boolDistance,'rows');
vUniqueNorm=v(:,iv);

% 
% flagKeepV=false(NV,1);
% for iV=1:NV
%     if find(boolDistance(iV,:),1,'first')==iV
%         flagKeepV(iV)=true;
%     end
% end
% vUniqueNorm=vNorm(:,flagKeepV);
