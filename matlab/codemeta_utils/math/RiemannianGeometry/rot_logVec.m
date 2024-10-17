function vVec=rot_logVec(R1,varargin)
if isempty(R1)
    R1=eye(3);
end
switch length(varargin)
    case 0
        R2=R1;
        R1=eye(3);
    otherwise
        if size(varargin{1},1)==3 && size(varargin{1},2)==3
            R2=varargin{1};
            varargin(1)=[];
        end
end
vVec=rot_vee(R1,rot_log(R1,R2,varargin{:}));
