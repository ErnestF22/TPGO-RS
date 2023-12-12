function result=localization_MLE_rigid_testSingle(t_node,varargin)
flagUseIsotropicInit=true;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flaguseisotropicinit'
            ivarargin=ivarargin+1;
            flagUseIsotropicInit=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

display('Without covariances')
t_nodeIsotropic=testNetworkAddDispersionMatricesRT(t_node,'methodR','identity','methodT','identity','methodCoupling','zero');
result.isotropic=runAndErros(t_nodeIsotropic);

display('With weights')
t_nodeWeighted=testNetworkAddDispersionMatricesRT(t_node,'givenIsotropic',t_node.dispersionMat);
if flagUseIsotropicInit
    t_nodeWeighted.gi=result.isotropic.t_node.gi;
end
result.weighted=runAndErros(t_nodeWeighted,'flagInit',flagUseIsotropicInit);

display('With covariances')
t_nodeCovariances=t_node;
if flagUseIsotropicInit
    t_nodeCovariances.gi=result.isotropic.t_node.gi;
end
result.covariances=runAndErros(t_nodeCovariances,'flagInit',flagUseIsotropicInit);

display('Spectral method')
result.spectral=runAndErros(t_node,'spectral');

function data=runAndErros(t_node,varargin)
flagSpectral=false;
flagInit=true;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flaginit'
            ivarargin=ivarargin+1;
            flagInit=varargin{ivarargin};
        case 'spectral'
            flagSpectral=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ~flagSpectral
    [t_node,errors]=localization_MLE_rigid(t_node,'flagInit',flagInit);
    data.intermediateErrors=errors;
else
    t_node=splitgij(t_node);
    Rij=t_node.Rij;
    Tij=t_node.Tij;
    E=t_node.E;
    t_node.Ri=sfm_rawAverageRotationsSpectral(Rij,E);
    t_node.Ti=sfm_rawAverageTranslationsDirect(t_node.Ri,Tij,E);
    t_node=mergegi(t_node);
end

t_node=testNetworkCompensate(t_node);
[rotErr,translErr]=testNetworkComputeErrors(t_node,'references');
rotErr=rotErr*180/pi;
translErr=translErr*180/pi;
disp(['Mean rot error    ' num2str(mean(rotErr))])
disp(['Std  rot error    ' num2str(std(rotErr))])
disp(['Mean transl error ' num2str(mean(translErr))])
disp(['Std  transl error ' num2str(std(translErr))])
%testNetworkShowErrors(t_node,'references')
data.t_node=t_node;
data.rotErr=rotErr;
data.translErr=translErr;
data.flagSpectral=flagSpectral;
data.flagInit=flagInit;
