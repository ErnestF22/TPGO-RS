function POCTagPositionError(TTruth,TOpt)
[d,TOptTransf]=procrustes(TTruth',TOpt','scaling',false,'reflection',false);
TOptTransf=TOptTransf';

disp('Procrustes error')
disp(d)
disp('TTruth/TOptTransf/Difference')
disp([TTruth;TOptTransf;TTruth-TOptTransf]);
