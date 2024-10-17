function m=locSegSegmentationPairIntersection(mA,mB)
m=zeros(size(mA));
%idx1=(mA==1 & mB>0) | (mB==1 & mA>0);
idx1=(mA==1 & mB==1);
m(idx1)=1;

