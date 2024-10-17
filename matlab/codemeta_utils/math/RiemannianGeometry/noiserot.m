%function R1=noiserot(R,sigmanoise);
%Given a rotation R, return a rotation R by perturbing the axis vee(R) with
%Gaussian noise with variance SIGMANOISE.
%If N is specified, NOISEROT creates an [3x3xN] array of noisy rotations
%If FLAGONEPARAMETER is set to true, only the norm of the axis is
%perturbed, resulting in a one parameter family of rotations

%%AUTORIGHTS%%

function R1=noiserot(R,sigmanoise,N,flagoneparameter,maxrot)
if(exist('N','var')==0)
    N=1;
end
if(exist('flagoneparameter','var')==0 || N==1)
    flagoneparameter=false;
end
if(exist('maxrot','var')==0)
    maxrot=pi/2;
end
if(size(R,3)==1)
    R1=zeros(3,3,N);
    switch (flagoneparameter)
        case 0
            for r=1:N
                Delta=cnormalize(randn(size(R,1),1));
                n=sigmanoise*randn;
                while(n>maxrot || n<-maxrot)
                    n=sigmanoise*randn;
                end           
                R1(:,:,r)=R*rot(Delta*n);
            end
        case 1
            x2=cnormalize(randn(size(R,1),1));
            for r=1:N
                t=sigmanoise*randn;
                while(t>maxrot)
                    t=sigmanoise*randn;
                end           
                R1(:,:,r)=R*rot(t*x2);
            end
        case 2
            x2=cnormalize(randn(size(R,1),1));
            x3=cnormalize(randn(size(R,1),1));
            for r=1:N
                t=sigmanoise*randn(2,1);
                while(max(t)>maxrot)
                    t=sigmanoise*randn(2,1);
                end           
                R1(:,:,r)=R*rot([x2 x3]*t);
            end
    end
else
    R1=zeros(size(R));
    for r=1:size(R1,3)
        R1(:,:,r)=noiserot(R(:,:,r),sigmanoise,1,flagoneparameter,maxrot);
    end
end
    
