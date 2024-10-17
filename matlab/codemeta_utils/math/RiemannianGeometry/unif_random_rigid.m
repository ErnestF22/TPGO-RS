%function Rarray=unif_random_rigid(N,mint,maxt)
%Compute an array of N uniformly distributed random rigid body
%transformations. The limits for the translation part are mint and maxt.
function Garray=unif_random_rigid(N,mint,maxt,maxtheta)
if(exist('maxtheta','var')==0)
    maxtheta=2*pi;
end
Garray=zeros(4,4,N);
for(it=1:N)
    Garray(1:3,1:3,it)=rot(cnormalize(randn(3,1))*maxtheta*(rand(1)));
    Garray(4,4,it)=1;
    Garray(1:3,4,it)=mint+rand(3,1).*(maxt-mint);
end