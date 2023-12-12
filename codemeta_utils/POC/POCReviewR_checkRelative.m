function POCReviewR_checkRelative
syms wij1 wij2 zij
syms wji1 wji2 zji

[~,~,Rij]=POCReviewRwz(wij1,wij2,zij);
%[~,~,Rji]=POCReviewRwz(wji1,wji2,zji);

[wji1b,wji2b]=POCReviewRzTow(Rij.');

disp(simplify(-wji1b))
disp(simplify(-wji2b))
