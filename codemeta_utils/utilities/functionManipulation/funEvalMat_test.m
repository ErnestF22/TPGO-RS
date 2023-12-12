function funEvalMat_test
X=randn(4,3,5);
f=@(x) x;
for k=1:length(size(X))
    disp(norm(vec(X-funEvalMat(f,X,'indexDimensionData',k))))
end