function R1=rot_retractions(R,A,type)
if exist('type','var')==0 || isempty(type)
    type='exp';
end

switch lower(type)
    case 'exp'
        R1=rot_exp(R,A);
    case 'qr'
        [R1,x]=qr(R+A);
%     case 'polar'
%         W=R'*A;
%         I=eye(size(R));
%         [U,S]=svd(I-W*W);
%         R1=R*(I+A)*U*diag(diag(S).^-0.5)*U';
end        
     