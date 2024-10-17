function [v,beta]=house(x)
%compute parameters for Householder transformation
%See Golub and Van Loan, p.210
sigma=x(2:end)'*x(2:end);
v=[1;x(2:end)];
if sigma==0
    %vector is in the direction of e1
    beta=0;
else
    %norm of x
    mu=sqrt(x(1)^2+sigma);
    if x(1)<=0
        v(1)=x(1)-mu;
    else
        v(1)=-sigma/(x(1)+mu);
    end
    beta=2*v(1)^2/(sigma+v(1)^2);
    v=v/v(1);
end