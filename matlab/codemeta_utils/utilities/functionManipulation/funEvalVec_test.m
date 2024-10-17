function funEvalVec_test
X=randn(4,3,5);
f=@(x) x;
for k=1:length(size(X))
    disp(norm(vec(X-funEvalVec(f,X,'indexDimensionVector',k))))
end