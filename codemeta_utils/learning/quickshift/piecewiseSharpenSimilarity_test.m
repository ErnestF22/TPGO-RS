function piecewiseSharpenSimilarity_test
A=[
    0.3860    1.0000         0         0    0.2105
    1.0000         0         0    0.2035    0.3805
         0         0    0.1964    0.8214    1.0000
    0.0090    0.1892    0.8198    1.0000         0
    0.1515    0.5758    1.0000         0    0.0303
];
threshold=median(A,2)-0.1;
AthresholdFull=piecewiseSharpenSimilarity(A,threshold);
AthresholdSparse=piecewiseSharpenSimilarity(sparse(A),threshold);

if ~issparse(AthresholdSparse)
    error('Result for sparse matrices should be sparse')
end
disp('Difference between full and sparse implementations')
disp(norm(AthresholdFull-AthresholdSparse,'fro'))

display(A)
display(threshold)
display(AthresholdFull)
