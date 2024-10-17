function rot3_expDiffSeries_test

skew=@(A) 0.5*(A-A');

A=skew(randn(3));

[expDiffSeriesDecomp(A) expDiffSeriesIt(A) rot3_expDiffSeries(A)]

function B=expDiffSeriesDecomp(A)
[V,S]=eig(A);
B=V*expDiffSeriesIt(S)*V';
;

function B=expDiffSeriesIt(A)
B=eye(size(A));
C=B;
for m=1:100
    C=C*(-A)/(m+1);
    B=B+C;
end

