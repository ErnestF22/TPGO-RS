function m=locSegSegmentationPairDiff(mA,mB)
m=zeros(size(mA));
idx1=mA==1 & (mB==0 | mB==2);
idx2=mA==2 & (mB==0 | mB==2);
m(idx1)=1;
m(idx2)=2;
