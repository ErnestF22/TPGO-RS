function acc15_presentation_example
load testData/testDataPresentation.mat
u=bearingCluster_getBearingsScalesFromE(x,E);
[m,L,s]=bearingCluster_clustering(E,u);

disp('L')
latex(L,'%.3f','nomath')

disp('L normalized')
latex(cnormalize(L')','%.3f','nomath')

