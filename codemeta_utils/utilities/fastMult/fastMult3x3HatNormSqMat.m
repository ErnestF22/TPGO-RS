%function B=fastMultHatNormSqMat(v,A)
%Same as fastMultMatHatNormSq, but with the order of the arguments reversed
%
%See also fastMultMatHatNormSq
function [B]=fastMult3x3HatNormSqMat(v,A)

c1=A(1:9:end).*v(1:3:end)+A(2:9:end).*v(2:3:end)+A(3:9:end).*v(3:3:end);
c2=A(4:9:end).*v(1:3:end)+A(5:9:end).*v(2:3:end)+A(6:9:end).*v(3:3:end);
c3=A(7:9:end).*v(1:3:end)+A(8:9:end).*v(2:3:end)+A(9:9:end).*v(3:3:end);

B=[...
    c1.*v(1:3:end);
    c1.*v(2:3:end);
    c1.*v(3:3:end);
    c2.*v(1:3:end);
    c2.*v(2:3:end);
    c2.*v(3:3:end);
    c3.*v(1:3:end);
    c3.*v(2:3:end);
    c3.*v(3:3:end);
   ];
B=B-reshape(A,9,[]);
