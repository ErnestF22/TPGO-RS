%Given a curve F(t), the function maps it to the tangent space at R using
%some basis and displays it
function rot3_tangentPlotFun(R,F,t,varargin)
style='b';
flagDisplayAxes=false;
flagShowVectorField=false;  %show a vector field (given as an opt parameter) along the curve
vfstep=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'axes'
            flagDisplayAxes=true;
        case 'vf'
            %Note: vf should give a tangent vector in the tangent space at F(t)
            ivarargin=ivarargin+1;
            vf=varargin{ivarargin};
            flagShowVectorField=true;
        case 'vfstep'
            ivarargin=ivarargin+1;
            vfstep=varargin{ivarargin};
        case 'style'
            ivarargin=ivarargin+1;
            style=varargin{ivarargin};
    end
    ivarargin=ivarargin+1;
end

hstate=ishold;
if flagDisplayAxes
    plot3([-pi pi],[0 0],[0 0],'k',[0 0],[-pi pi],[0 0],'k',[0 0],[0 0],[-pi pi],'k',0,0,0,'o');
    hold on
end

f=@(t) rot_vee(R,rot_log(R,F(t)));
[ft,t]=plotfun(f,t,style,'plotDim',3);

if flagShowVectorField
    tvf=t(1:vfstep:end);
    Nt=length(tvf);
    v=zeros(3,Nt);
    for it=1:Nt
        Ft=F(tvf(it));
        v(:,it)=rot_vee(R,rot_parallel(Ft,rot_log(Ft,R),vf(tvf(it))));
    end
    quiver3(ft(1,1:vfstep:end),ft(2,1:vfstep:end),ft(3,1:vfstep:end),v(1,:),v(2,:),v(3,:))
end

if flagDisplayAxes && ~hstate
    hold off
end
