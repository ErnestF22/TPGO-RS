%function rot3_plotvf(S,R,vf,varargin)
%Displays the vector field vf at points R by mapping all in the tangent
%space at S

%With "vee representation" we intend the operation of passing from rotation
%to points in the tangent space at a origin rotation S to vectors obtained
%as coefficients by fixing a basis for the tangent space at S
function rot3_plotvf(S,R,vf,varargin)
flagVeeRepresentation=false;    %both R and vf are already given in vector form (no need to take the vee)
flagDisplayAxes=false;
flagDisplayPoints=false;
flagTangentMapping='parallel';   %determines how to transport vectors between tangent spaces: parallel or pushforward
flagStats=false;

I1=[]; I2=[]; I3=[];
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'axes'
            flagDisplayAxes=true;
        case 'points'
            flagDisplayPoints=true;
        case 'vee'
            flagVeeRepresentation=true;
        case 'mapping'
            ivarargin=ivarargin+1;
            flagTangentMapping=varargin{ivarargin};
        case 'intervals'
            ivarargin=ivarargin+1;
            I1=varargin{ivarargin};
            ivarargin=ivarargin+1;
            I2=varargin{ivarargin};
            ivarargin=ivarargin+1;
            I3=varargin{ivarargin};
        case 'stats'
            flagStats=true;
    end
    ivarargin=ivarargin+1;
end

%automatic generation of R on a "grid", if requested
while ischar(R)
    switch R
        case 'grid'
            R='7x7x11';
        otherwise
            [N,nN]=sscanf(R,'%dx%dx%d');
            if nN<3
                error('The string for automatically generating a grid must be in the format ''N1xN2xN3''')
            end
            R=genRotGrid(N,I1,I2,I3);
    end
end
  

%get number of locations at which we need to show vf
if ~flagVeeRepresentation
    Np=size(R,3);
else
    Np=size(R,2);
end

%if vf is a function, evaluate it at all the locations
if isa(vf,'function_handle')
    if flagVeeRepresentation
        error('Cannot evaluate vf with vee representation (not implemented yet!)')
    end
    v=zeros(3,3,Np);
    for it=1:Np
        v(:,:,it)=vf(R(:,:,it));
    end
    vf=v;
    clear v;
end

%pass to vee representation
if ~flagVeeRepresentation
    Np=size(vf,3);
    for it=1:Np
        R1=R(:,:,it);
        switch flagTangentMapping
            case 'parallel'
                vf(:,:,it)=rot_parallel(R1,rot_log(R1,S),vf(:,:,it));
            case 'pushforward'
                vf(:,:,it)=S*R1'*vf(:,:,it);
            otherwise
                error('Type of mapping not recognized')
        end
    end
    vf=rot_vee(S,vf);
    R=rot_vee(S,rot_log(S,R));
end

%display
hstate=ishold;
if flagDisplayAxes
    plot3([-pi pi],[0 0],[0 0],'k',[0 0],[-pi pi],[0 0],'k',[0 0],[0 0],[-pi pi],'k',0,0,0,'o');
    hold on
end

if flagDisplayPoints
    plot3(R(1,:),R(2,:),R(3,:),'.');
    hold on
end

quiver3(R(1,:),R(2,:),R(3,:),vf(1,:),vf(2,:),vf(3,:))
axis equal
set(gcf,'renderer','opengl')

if (flagDisplayAxes || flagDisplayPoints) && ~hstate
    hold off
end

if flagStats
    normvf=rot_tangentNorm(S,rot_hat(S,vf));
    disp(['v.f. norm (MIN/MAX): ' num2str(min(normvf)) '/'  num2str(max(normvf))]);
end

function R=genRotGrid(N,I1,I2,I3)
%N(1), alpha, I1 - distance base rotation from origin rotation
%N(2), beta, I2 - angle in the plane perpendicular to direction to the origin
%N(3), gamma, I3 - number of points on the geodesic
epsilon=1e-6;

if ~exist('I1','var') || isempty(I1)
    I1=linspace(-pi+epsilon,pi-epsilon,N(1));
end

if ~exist('I2','var') || isempty(I2)
    I2=linspace(0,pi,N(2));
end

if ~exist('I3','var') || isempty(I3)
    I3=linspace(-pi+epsilon,pi-epsilon,N(3));
end


R=[];
v0=hat([1;0;0]);
for alpha=I1
    R0=rot_exp(eye(3),alpha*v0);
    
    for beta=I2
        vRId=cos(beta)*hat([0;1;0])+sin(beta)*hat([0;0;1]);
        
        v=rot_parallel(eye(3),alpha*v0,vRId);        
        Rt=@(t) rot_exp(R0,t*v);
        for gamma=I3
            R=cat(3,R,Rt(gamma));
        end
    end
end
