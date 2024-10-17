%Compute the householder rotation that maps one vector direction to another
%function H=householderRotation(x1,x2)
%This function has two modes of operation, depending on the size of x2
%Mode 1: size(x2,1)==1
%H will map x1 to the x2-th standard vector.
%This is a modified version of algorithm 5.1.1 at on page 210 of Golub and
%Van Loan, 3rd ed. The algorithm below first flips the sign of x and then
%flips the sign of H. In this way, the final transformation H will be a
%rotation instead of a reflection. Also, it permits to specify the final
%direction to be any standard basis vector ek instead of e1.
%Mode 2: size(x2,1)==size(x1,1)
%H will map an arbitrary x1 to another arbitrary x2. This will use the
%algorithm in the notes. Note that, in principle, this is numerically less
%accurate than the previous mode when x1 and x2 are nearly opposite
function H=householderRotation(x1,x2)
if ~exist('x2','var')
    x2=1;
end

[D1]=size(x1,1);
[D2]=size(x2,1);

[N,x1,x2]=argRepmat(2,x1,x2);

if N>1
    H=zeros(D1,D1,N);
    for iN=1:N
        H(:,:,iN)=householderRotation(x1(:,iN),x2(:,iN));
    end
else
    flagStandardVector=(D2==1);

    if flagStandardVector
        x1=-x1;

        n=length(x1);
        xRed=x1([1:x2-1 x2+1:n]);
        sigma=xRed'*xRed;
        v=x1;
        v(x2)=1;

        if sigma==0
            beta=0;
        else
            mu=sqrt(x1(x2)^2+sigma);
            if x1(x2)<=0
                v(x2)=x1(x2)-mu;
            else
                v(x2)=-sigma/(x1(x2)+mu);
            end
            beta=2*v(x2)^2/(sigma+v(x2)^2);
            v=v/v(x2);
        end
        H=-eye(D1)+beta*(v*v');
    else
        v=cnormalize(cnormalize(x1)+cnormalize(x2));
        H=-eye(D1)+2*(v*v');
    end
end
