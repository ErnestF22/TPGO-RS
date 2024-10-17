function Y1=grassman_retractions(Y,H,type)
if exist('type','var')==0 || isempty(type)
    type='exp';
end

switch lower(type)
    case 'exp'
        Y1=grassman_exp(Y,H);
    case 'qr'
        [Y1,void]=qr(Y+H);
        Y1=Y1(:,1:size(Y,2));
    otherwise
        error('Retraction type unknown')
end        
     