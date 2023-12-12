%function R=prodrot(Rarray)
%Computes the product Rarray(:,:,N)*Rarray(:,:,N-1)*...*Rarray(:,:,1)
function R=prodrot(Rarray,opt)
if(exist('opt','var')==0)
    opt='dir';
end
R=Rarray(:,:,1);
switch(opt)
    case 'dir'
        for(r=2:size(Rarray,3))
            R=Rarray(:,:,r)*R;
        end
    case 'rev'
        for(r=2:size(Rarray,3))
            R=R*Rarray(:,:,r);
        end
end
