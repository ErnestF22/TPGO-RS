function POCNumericalDerivatives
filterOrder=2;
filterWindow=41;
switch 2
    case 1
        load('sampleDynamicSfMDataset_interpolation')

        [xbFilter,dxbFilter,ddxbFilter,g]=sgolayFilterDerivatives(mean(diff(t)),xb,filterOrder,filterWindow);
        idxValid=filterWindow:length(t)-filterWindow;
        dxbConv=convVectorized(dxb,g(:,1),'same');
        subsample idxValid xb xbFilter dxb dxbFilter ddxb ddxbFilter ...
            Gammab wIMU alphaIMU Rbs taub nub dxbConv

        disp('RMSE xbFilter dxbFilter ddxbFilter')
        disp([rmse(xbFilter, xb) rmse(dxbFilter, dxb) rmse(ddxbFilter, ddxb)])
        disp('RMSE dxbFilter-dxb ddxbFilter-dxbConv')
        disp([rmse(dxbFilter, dxb) rmse(dxbFilter, dxbConv)])
    case 2
        %flagHighFreqDisturbance=false;
        flagHighFreqDisturbance=true;
        t=linspace(1,2,400);
        TSampling=mean(diff(t));
        s=sin(2*pi*t);
        ds=2*pi*cos(2*pi*t);
        if flagHighFreqDisturbance
            s=s+sin(40*pi*t);
            ds=ds+40*pi*cos(40*pi*t);
        end
        funCheckDerInterpInterp(t,s,ds,t)
        [~,g]=sgolay(filterOrder,filterWindow);
        dsEst=-conv(s,g(:,2),'same')/TSampling;
        dsFilt=conv(ds,g(:,1),'same');
        idxValid=filterWindow:length(t)-filterWindow;
        plot(t(idxValid),dsEst(idxValid),t(idxValid),dsFilt(idxValid),':')
end

function subsample(nameIdxVar,varargin)
for ivarargin=1:length(varargin)
    var=varargin{ivarargin};
    A=evalin('caller',var);
    sz=size(A);
    if sz(2)==1
        sz(2)=[];
    end
    switch length(sz)
        case 1
            evalin('caller',[var '=' var '(' nameIdxVar ');']);
        case 2
            evalin('caller',[var '=' var '(:,' nameIdxVar ');']);
        case 3
            evalin('caller',[var '=' var '(:,:,' nameIdxVar ');']);
    end
end
