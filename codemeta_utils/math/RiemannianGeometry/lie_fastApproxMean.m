function ymean=lie_fastApproxMean(y,lie_funs)
    N=size(y,3);
    ymean=y(:,:,1);
    for n=2:N
        ymean=lie_funs.exp(ymean,(n-1)/n*lie_funs.log(ymean,y(:,:,n)));
    end
end
