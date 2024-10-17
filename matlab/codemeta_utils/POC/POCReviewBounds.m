function POCReviewBounds
switch 3
    case 1
        f=@(x,eta) [
            0.5+0.5*x;
            0.5+(0.5-eta)*x;
            (1-eta)*x;
            ];
        funPlot(@(eta) f(0.25,eta));
        legend('f0','f1','fb')

    case 2
        syms delta w epsilon d eta
        
        %boundP=@(x) log((1-eta)*exp(-x^2));
        boundP=@(x) log((1-eta))-x^2;
        %boundM=@(x) log(eta*exp(-x^2));
        boundM=@(x) log(eta)-x^2;
        multP=1/2*(1+w);
        multM=1/2*(1-w);
        Hij1=-multP*boundP(d)-multM*boundM(d);
        Hij0=-multP*boundP(epsilon)-multM*boundM(epsilon);
        Hij=simplify(delta*w*Hij1+(1-delta*w)*Hij0);
        HijPaper=delta*w*d^2+(1-delta*w)*epsilon^2+w/2*(log(eta)-log(1-eta))-1/2*(log(eta)+log(1-eta));
        keyboard
        
    case 3
        syms delta w w0 epsilon d eta lambda
        J=lambda*(delta*w*d^2+(1-delta*w)*epsilon^2+w*eta)+(1-lambda)*(w-w0)^2;
        solve(diff(J,w)==0,w)
end

       


