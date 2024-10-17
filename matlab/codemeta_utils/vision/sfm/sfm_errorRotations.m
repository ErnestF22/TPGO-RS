function e=sfm_errorRotations(data,varargin)
memberAbsolutePoses='poseTruth';
memberRelativePoses='matchPoseTruth';
memberMatch='matchFiltered';

methodAbsolutePoses='reference';
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'memberabsoluteposes'
            ivarargin=ivarargin+1;
            memberAbsolutePoses=varargin{ivarargin};
        case 'memberrelativeposes'
            ivarargin=ivarargin+1;
            memberRelativePoses=varargin{ivarargin};
        case 'membermatch'
            ivarargin=ivarargin+1;
            memberMatch=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

Ri=G2R(data.(memberAbsolutePoses));
Rij=G2R(data.(memberRelativePoses));

E=(sfm_getMatchIdxImg(data,memberMatch))';
NE=size(E,1);
e=zeros(NE,1);
for iE=1:NE
    iCamera=E(iE,1);
    jCamera=E(iE,2);
    e(iE)=rot_dist(Rij(:,:,iE), Ri(:,:,iCamera)'*Ri(:,:,jCamera));
end

if nargout==0
    plot(e*180/pi)
    xlabel('Edge index')
    ylabel('Rotation error [deg]')
end
