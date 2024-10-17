%function draw3dcameraFromRT(R,T,varargin)
%Draw a pyramid and a set of axes to represent a camera
%Optional argument
%   'FlipZ'         Flip the z axis when drawing the pyramid
%   'Color1',c      Color for the sides of the pyramid
%   'Color2',c      Color for the base of the pyramid
%   'Scale',s       Scale for the pyramid (if 'Scale' is not specified, s=1)
%   'flagAxes',f    Flag to specify if axes should be drawn or not
%   'flagAlpha',f   Flag to enable/disable reduced opacity
%   'references'    R,T are in the "reference" interpretation
%   'poses'         R,T are in the "pose" interpretation
%
%See also draw3dcameraFromPose

%%AUTORIGHTS%%
function draw3dcameraFromRT(R,T,varargin)

if(nargin<1)
    R=eye(3);
end
if(nargin<2)
    T=[0;0;0];
end

NPoses=size(R,3);
if NPoses>1
    flagHold=ishold();
    for iPose=1:NPoses
        draw3dcameraFromRT(R(:,:,iPose),T(:,iPose),varargin{:});
        hold on
    end
    if ~flagHold
        hold off
    end
else
    methodAbsolutePoses='reference';
    shape='camera';

    flagAlpha=true;
    flagFlipZ=false;
    flagAxes=true;
    flagAxesLabels=true;
    flagPlotVertexNumber=false;
    color1=[15934	35723	14392]/65535/0.6;
    color2=[514	13107	39321]/65535;
    scale=1;
    alpha=0.7;

    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'posttranslation'
                flagpostT=true;
            case 'color1'
                ivarargin=ivarargin+1;
               color1=varargin{ivarargin};
            case 'color2'
                ivarargin=ivarargin+1;
               color2=varargin{ivarargin};
            case 'color'
                ivarargin=ivarargin+1;
               color1=varargin{ivarargin};
               color2=color1;
            case 'alpha'
                ivarargin=ivarargin+1;
                alpha=varargin{ivarargin};
            case 'alpha'
                ivarargin=ivarargin+1;
                alpha=varargin{ivarargin};
            case 'flagalpha'
                ivarargin=ivarargin+1;
                flagAlpha=varargin{ivarargin};
            case 'flipz'
                flagFlipZ=true;
            case 'scale'
                ivarargin=ivarargin+1;
                scale=varargin{ivarargin};
            case 'flagaxes'
                ivarargin=ivarargin+1;
                flagAxes=varargin{ivarargin};
            case 'flagaxeslabel'
                ivarargin=ivarargin+1;
                flagAxesLabels=varargin{ivarargin};
            case 'nolabels'
                flagAxesLabels=false;
            case 'references'
                methodAbsolutePoses='reference';
            case 'poses'
                methodAbsolutePoses='pose';
            case 'methodabsoluteposes'
                ivarargin=ivarargin+1;
                methodAbsolutePoses=lower(varargin{ivarargin});
            case 'shape'
                ivarargin=ivarargin+1;
                shape=lower(varargin{ivarargin});
            otherwise
                error('Argument not valid!')
        end
        ivarargin=ivarargin+1;
    end

    switch length(scale)
        case 1
            scale=ones(3,1)*scale;
        case 2
            scale(3)=1;
        case 3
            %nothing to do
        otherwise
            error('Incorrect scale dimensions')
    end

    flagShowPlane=false;
    switch shape
        case 'camera'
            x=[ 0  1  1 -1 -1];
            y=[ 0  1 -1 -1  1];
            z=[ 0  1  1  1  1];

            F=[  1 2 3
                 1 3 4;...
                 1 4 5;...
                 1 5 2]';
            flagShowPlane=true;
        case 'velodyne'
            x=0.5*[ 0  1  1 -1 -1  0  1  1 -1 -1];
            y=0.5*[ 0  1 -1 -1  1  0  1 -1 -1  1];
            z=[ 1  1  1  1  1 -1 -1 -1 -1 -1];

            F=[  1 2 3;
                 1 2 5;
                 1 4 5;
                 6 7 8;
                 6 7 10;
                 6 9 10;
                 6 3 8;
                 6 3 1;
                 6 4 1;
                 6 4 9;
                 5 7 10;
                 5 7 2]';
        case 'hokuyo'
            x=[ 0  0  0];
            y=[ 0  1  1];
            z=[ 0 -1  1];

            F=[  1 2 3]';
        case 'plane'
            %no shape defined
        case 'april'
            x=2*([0 1 1 0]-0.5);
            y=2*([0 0 1 1]-0.5);
            z=[0 0 0 0];
            F=[ 1 2 3;
                1 4 3]';
        otherwise
            error('Shape name not recognized')
    end
    if flagFlipZ
        z=-z;
    end

    flaghold=ishold;

    if strcmpi(methodAbsolutePoses,'pose')
        [R,T]=invRT(R,T);
    end

    optsPlotFrame={'scale',scale};

    if flagAxes
        if flagAxesLabels
            plot3dframe(T,R,[],optsPlotFrame{:})
        else
            plot3dframe(T,R,[],'nolabels',optsPlotFrame{:})
        end
    end
    hold on
    if ~strcmpi(shape,'plane')
        %process vertices information
        V=diag(scale)*[x;y;z];

        if flagPlotVertexNumber
            for iV=1:size(V,2);
                text(V(1,iV),V(2,iV),V(3,iV),num2str(iV));
            end
            hold on
        end

        V=R*V+T*ones(1,size(V,2));

        %draw polygon
        h=draw3dpolygon(V,F);
        if flagAlpha
            set(h,'FaceAlpha',alpha)
        end
        set(h,'FaceColor',color1)
    end
    if flagShowPlane
        h=patch(V(1,2:end),V(2,2:end),V(3,2:end), 'b');
        if flagAlpha
            set(h,'FaceAlpha',0.7)
        end
        set(h,'FaceColor',color2)
    end

    if(~flaghold)
        hold off
    end
end