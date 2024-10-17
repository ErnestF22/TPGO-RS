function bearingCluster_tikz_test
fid=fopen('tikzTest/plot.tex','wt');

[A,x]=bearingCluster_generateTest('butterfly-bend');
E=adj2edges(A,'oriented');

bearingCluster_tikz(x,E,'fileId',fid);

fclose(fid);

cdOld=cd('tikzTest');
system('pdflatex head.tex');
cd(cdOld);
