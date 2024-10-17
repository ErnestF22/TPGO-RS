function clipToBox_test
x=polarGrid(1,21);
xClip=clipToBox(x,[4.5 4]);
plotEdgesOneToMany(zeros(2,1),x,'k')
hold on
plotPoints(xClip,'rx')
hold off

axis equal
