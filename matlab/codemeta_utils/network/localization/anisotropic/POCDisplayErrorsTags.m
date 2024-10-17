function POCDisplayErrorsTags(t_node,color)
if nargin<2
    color='r';
end
TOpt=subMatrix(G2T(cat(3,t_node.gi)),1:3,1:6);
TTruth=subMatrix(G2T(cat(3,t_node.gitruth)),1:3,1:6);
[d,TOptTransf]=procrustes(TTruth',TOpt','scaling',false,'reflection',false);
TOptTransf=TOptTransf';

disp('Sqrt procrustes error')
disp(sqrt(d))
% disp('TTruth/TOptTransf/Difference')
% disp([TTruth;TOptTransf;TTruth-TOptTransf]);

plotPoints(TTruth)
hold on
plotPoints(TOptTransf,{'Color',color})
hold off
axis equal
