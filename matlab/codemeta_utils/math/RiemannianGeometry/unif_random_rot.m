%function Rarray=unif_random_rot(N)
%Compute an array of N uniformly distributed random rotations

%%AUTORIGHTS%%

function Rarray=unif_random_rot(N,varargin)
B=eye(3);
maxangle=2*pi;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'subfamily'
            ivarargin=ivarargin+1;
            B=orth(varargin{ivarargin});
            B=B*B';
        case 'maxangle'
            ivarargin=ivarargin+1;
            maxangle=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

Rarray=zeros(3,3,N);
for it=1:N
    Rarray(:,:,it)=rot(B*cnormalize(randn(3,1))*maxangle*(rand(1)));
end