%function B=fastMultMatHat(A,v)
%Given a [3x3xN] array A and a [3xN] matrix v, computes the array B such
%that B(:,n)=reshape(A(:,:,n)*hat(v(:,n)),[],1)
%The computations are vectorized along the third dimension and make use of
%the fact that some entries of hat(v) are zero. This function is especially
%faster than the naive implementation with a for loop for large N (N>50)
%If you want B in the form of a [3x3xN] array, use 
% B=reshape(B,3,3,[])
function B=fastMult3x3MatHat(A,v)

B=[...
    -A(4:9:end).*v(3:3:end)+A(7:9:end).*v(2:3:end);
    -A(5:9:end).*v(3:3:end)+A(8:9:end).*v(2:3:end);
    -A(6:9:end).*v(3:3:end)+A(9:9:end).*v(2:3:end);
    A(1:9:end).*v(3:3:end)-A(7:9:end).*v(1:3:end);
    A(2:9:end).*v(3:3:end)-A(8:9:end).*v(1:3:end);
    A(3:9:end).*v(3:3:end)-A(9:9:end).*v(1:3:end);
    -A(1:9:end).*v(2:3:end)+A(4:9:end).*v(1:3:end);
    -A(2:9:end).*v(2:3:end)+A(5:9:end).*v(1:3:end);
    -A(3:9:end).*v(2:3:end)+A(6:9:end).*v(1:3:end);
    ];
