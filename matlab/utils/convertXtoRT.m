function X = convertXtoRT(Xvec, p, d, n)
    XRvec = Xvec(1:p*d*n);
    XTvec = Xvec(p*d*n + 1:end);
    XRhst = reshape(XRvec, p, []);
    X.R = matUnstackH(XRhst, d);
    X.T = reshape(XTvec, p, n);
end