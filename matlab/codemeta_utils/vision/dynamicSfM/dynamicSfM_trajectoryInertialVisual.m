function dynamicSfM_trajectoryInertialVisual(fileName)
fileNameLoad=fileName;
fprintf('Loading file to %s\n',fileNameLoad);
load(fileNameLoad)

%methodDecimation='original';
%methodDecimation='sampling';
methodDecimation='interpolation';

switch methodDecimation
    case 'original'
        %Nothing to do
    case 'sampling'
        idxSampling=[1:10:length(t)];
        subsample Rsb Rbs Tsb dTsb ddTsb wb taub nub dwb alphab wIMU alphaIMU Gammab Theta t
    case 'interpolation'
        tQuery=min(t):1/30:max(t);
        Rsb=rot_interpolationLinear(t,Rsb,tQuery);
        Rbs=rot_interpolationLinear(t,Rbs,tQuery);
        interpolate Tsb dTsb ddTsb wb taub nub dwb alphab wIMU alphaIMU Gammab Theta
        t=tQuery;
    otherwise
        error('Method decimation not recognized');
end

SX=4;
SY=3;
SZ=2;
X=[repmat(reshape(permute(cat(3,repmat(1:SX,SY,1),repmat((1:SY)',1,SX)),[3 1 2]),2,[]),1,SZ); kron(1:SZ,ones(1,SX*SY))];
X=center(X);

%Coordinates of the geometry (3-D points) in the camera frame and their
%derivatives
Xb=rigidTransform(Rsb,Tsb,X,'references','wc');
dXb=dynSfM_derGeometry(Xb,wb,nub);
ddXb=dynSfM_derDerGeometry(Xb,wb,dwb,nub,alphab);

%Ground truth
XbTruth=Xb;
taubTruth=taub;
nubTruth=nub;
RbsTruth=Rbs;

%Image measurements
P=[eye(2) zeros(2,1)];
xb=multiprod(P,Xb);
dxb=multiprod(P,dXb);
ddxb=multiprod(P,ddXb);

%Ground truth factorization
MTruth=Rbs;
STruth=X;
TTruth=permute(taub,[1 3 2]);
TpTruth=permute(nub,[1 3 2]);
TppTruth=permute(alphab,[1 3 2]);
disp(' - Consistency check on factorization')
NX=size(X,2);
disp(norm(vec(matStack(XbTruth)-matStack(MTruth)*matStack(STruth)+matStack(multiprod(TTruth,ones(1,NX)))),'inf'))


fileNameSave=[fileName 'Visual_' methodDecimation];
save(fileNameSave,...
    'XbTruth','taubTruth','nubTruth','RbsTruth',... %Ground truth reconstruction
    'xb','dxb','ddxb','alphaIMU','wIMU',... %Image and IMU measurements
    'MTruth','STruth','TTruth','TpTruth','TppTruth',... %Ground truth factorization
    'Rbs','Tsb','dTsb','ddTsb','X','taub','nub',... %Ground truth trajectory
    'Gammab','Theta',... %Known inputs
    'g','eThrust','J','m','P','Rbc',... %known parameters of the system
    't'... %Sampling times
    )
fprintf('Dataset saved to %s\n',fileNameSave);

figure(1)
plotPoints(Tsb,'-')
hold on
plotPoints(X,'x')
draw3dcameraFromRT(Rsb(:,:,1),Tsb(:,1),'scale',0.2)
hold off
view(0,90)
axis equal

function subsample(varargin)
for ivarargin=1:length(varargin)
    var=varargin{ivarargin};
    A=evalin('caller',var);
    sz=size(A);
    if sz(2)==1
        sz(2)=[];
    end
    switch length(sz)
        case 1
            evalin('caller',[var '=' var '(idxSampling);']);
        case 2
            evalin('caller',[var '=' var '(:,idxSampling);']);
        case 3
            evalin('caller',[var '=' var '(:,:,idxSampling);']);
    end
end

function interpolate(varargin)
for ivarargin=1:length(varargin)
    var=varargin{ivarargin};
    evalin('caller',[var '=real_interp(t,' var ',tQuery,''spline'');']);
end