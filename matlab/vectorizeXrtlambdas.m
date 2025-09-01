function Xvec = vectorizeXrtlambdas(X)
    Xvec = [X.R(:); X.T(:); X.lambda(:)];
    
end