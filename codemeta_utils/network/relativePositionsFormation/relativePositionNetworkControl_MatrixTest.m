function relativePositionNetworkControl_MatrixTest
resetRands();
d=2;

NNodes=5;
ECost=[...
     1     2
     2     1
     1     3
     3     1
     2     3
     3     2
     2     4
     4     2
     3     4
     4     3
     3     5
     5     3
     4     5
     5     4
     ];
t_node=relativePositionNetworkTestNetwork(NNodes);
t_node=relativePositionNetworkTestNetwork(NNodes);
TiTruth=t_node.TiTruth;
%EControl=ECost([1:10 12 14],:);
EControl=[ECost(ECost(:,1)>ECost(:,2),:); 1 3; 3 5];
TijTruth=relativePositionNetworkCompute(TiTruth,ECost);
TijControlTruth=relativePositionNetworkCompute(TiTruth,EControl);

Ti=5*randn(d,NNodes);
Tij=relativePositionNetworkCompute(Ti,ECost);
TijControl=relativePositionNetworkCompute(Ti,EControl);

dx=relativePositionNetworkControl(EControl,TijControl,TijControlTruth);
dxEdges=relativePositionNetworkControlEdges(EControl,TijControl,TijControlTruth);

B=edges2incmatrix(fliplr(EControl),NNodes,'directed');
Bd=kron(B,eye(d));
disp(TijControl-reshape(Bd*Ti(:),size(TijControl)))
disp(TijControlTruth-reshape(Bd*TiTruth(:),size(TijControl)))
disp(dxEdges-reshape([Bd -Bd]*[Ti(:);TiTruth(:)],size(TijControl)))
disp(dx-reshape((Bd<0)'*dxEdges(:),d,NNodes))
A=(B<0)'*B;
Ad=kron(A,eye(2));
disp(dx-reshape([Ad -Ad]*[Ti(:);TiTruth(:)],d,NNodes))
