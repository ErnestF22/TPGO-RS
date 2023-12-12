%Computes the inner product between packed tangent vectors
%function c=rotDyn_metricPacked(v1,v2)
function c=rotDyn_metricPacked(v1,v2)
v1(1:9,:)=v1(1:9,:)/2;
c=squeeze(multiprod(multitransp(permute(v1(1:12,:),[1 3 2])),permute(v2(1:12,:),[1 3 2])))';
