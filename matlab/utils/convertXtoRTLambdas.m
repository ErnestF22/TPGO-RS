function X = convertXtoRTLambdas(Xvec, p, d, n)
    
    XRvec = Xvec(1:p*d*n);
    XTvec = Xvec(p*d*n + 1:p*d*n + p*n);
    XRhst = reshape(XRvec, p, d*n);
    X.R = matUnstackH(XRhst, d);
    X.T = reshape(XTvec, p, n);
    X.lambda = Xvec(p*d*n + p*n + 1:end);

end

