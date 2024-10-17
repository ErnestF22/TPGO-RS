%function Y1=stiefel_retractions(Y,H,type)
%Compute a retraction for the Stiefel manifold from Y in the direction H
%The argument type can be 'polar' (default), 'qr', or 'exp'
function Y1=stiefel_retractions(Y,H,type)
if exist('type','var')==0 || isempty(type)
    type='polar';
end

switch lower(type)
    case 'exp'
        Y1=stiefel_exp(Y,H);
    case 'qr'
        [Y1,void]=qr(Y+H);
        Y1=Y1(:,1:size(Y,2));
    case 'polar'
        [V,S]=eig(eye(size(Y,2))+H'*H);
        Y1=(Y+H)/(V*diag(sqrt(diag(S)))*V');
        
    otherwise
        error('Retraction type unknown')
end        
     