function m=stiefel_metric(Y,H1,H2,method)
if ~exist('method','var')
    method='euclidean';
end



[n,p]=size(Y);
Inp=[eye(p); zeros(n-p,p)];
switch lower(method)
    case 'canonical'
        %m=trace(H1'*(eye(n)-0.5*Y*Y')*H2);
        P=repmat(eye(n),[1,1,size(Y,3)])-0.5*multiprod(Y,multitransp(Y));
        m=multitrace(multiprod(multitransp(H1),multiprod(P,H2)));
    case 'euclidean'
        m=multitrace(multiprod(multitransp(H1),H2));
    otherwise
        error('Invalid method for Stiefel metric')
end
end