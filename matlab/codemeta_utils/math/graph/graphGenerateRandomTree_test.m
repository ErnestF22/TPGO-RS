function graphGenerateRandomTree_test
%generate a complete graph
V=2:7;
[R,C]=meshgrid(V,V);
E=[R(:),C(:)];

[ETree,VTree]=graphGenerateRandomTree(E);

display(VTree)
display(ETree)
