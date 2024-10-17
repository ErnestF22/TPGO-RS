function th=sfm_rawThresholdEstimate(e)
w=2;
m=median(e);
md=median(abs(e-m));
th=m+w*md;
