%function C=fastMult3x3MatMat(A,B,C)
%Given two [3x3xN] arrays A and B, computes the matrix C such
%that C(:,n)=reshape(A(:,:,n)*B(:,:,n),[],1). The matrices A,B can also be
%reshaped (e.g., vectorized).
%The computations are vectorized along the third dimension. This function
%is significantly faster than the naive implementation with a for loop for
%large N (N>75).
%If you want C in the form of a [3x3xN] array, use 
%   C=reshape(C,3,3,[])
%
%Optionally, one can pass in a placeholder for the output arguments to
%avoid new memory allocations (tested on R2009b)

%%AUTORIGHTS%%
function C=fastMult3x3MatMat(A,B,C)

L=numel(A);

C=[...
    A(1:9:L).*B(1:9:L)+A(4:9:L).*B(2:9:L)+A(7:9:L).*B(3:9:L);
    A(2:9:L).*B(1:9:L)+A(5:9:L).*B(2:9:L)+A(8:9:L).*B(3:9:L);
    A(3:9:L).*B(1:9:L)+A(6:9:L).*B(2:9:L)+A(9:9:L).*B(3:9:L);
    A(1:9:L).*B(4:9:L)+A(4:9:L).*B(5:9:L)+A(7:9:L).*B(6:9:L);
    A(2:9:L).*B(4:9:L)+A(5:9:L).*B(5:9:L)+A(8:9:L).*B(6:9:L);
    A(3:9:L).*B(4:9:L)+A(6:9:L).*B(5:9:L)+A(9:9:L).*B(6:9:L);
    A(1:9:L).*B(7:9:L)+A(4:9:L).*B(8:9:L)+A(7:9:L).*B(9:9:L);
    A(2:9:L).*B(7:9:L)+A(5:9:L).*B(8:9:L)+A(8:9:L).*B(9:9:L);
    A(3:9:L).*B(7:9:L)+A(6:9:L).*B(8:9:L)+A(9:9:L).*B(9:9:L);
   ];
