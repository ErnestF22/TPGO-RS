%Compute the residuals of points from an essential matrix
function e=epipolarResiduals(E,x1,x2,type)
if ~exist('type','var') || isempty(type)
    type='algebraicSq';
end
NPoints=size(x1,2);
NE=size(E,3);

e=zeros(NE,NPoints);
for iE=1:NE
    Ei=E(:,:,iE);
    switch lower(type)
        case 'algebraicabs'
            ei=abs(epipolarConstraintFromE(Ei,x1,x2));
        case 'algebraicsq'
            ei=epipolarConstraintFromE(Ei,x1,x2).^2;
        case 'linedistanceabs'
            ei=sqrt(epipolarResidualsLineDistanceSq(Ei,x1,x2));
        case 'linedistancesq'
            ei=epipolarResidualsLineDistanceSq(Ei,x1,x2).^2;
        case 'sampsonabs'
            ei=sqrt(epipolarResidualsSampsonSq(Ei,x1,x2));
        case 'sampsonsq'
            ei=epipolarResidualsSampsonSq(Ei,x1,x2);
        otherwise
            error('Type of epipolar residuals not recognized');
    end
    e(iE,:)=ei;
end
