function [Q1,v1]=essential_exp(Q,v,varargin)

flagProject=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
%         case 'method'
%             ivarargin=ivarargin+1;
%             method=varargin{ivarargin};
        case 'project'
            flagProject=true;

        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagProject
    v=essential_tangentProj(Q,v);
end

Q1(1:3,:,:)=rot_exp(Q(1:3,:,:),v(1:3,:,:));
Q1(4:6,:,:)=rot_exp(Q(4:6,:,:),v(4:6,:,:));
if nargout>1
    v1(1:3,:,:)=multiprod(Q1(1:3,:,:),multiprod(multitransp(Q(1:3,:,:)),v(1:3,:,:)));
    v1(4:6,:,:)=multiprod(Q1(4:6,:,:),multiprod(multitransp(Q(4:6,:,:)),v(4:6,:,:)));
end