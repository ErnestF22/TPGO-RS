function POCReviewR_check24
wij=sym('wij',[2 1]);
twij=[wij;0];
omegai=sym('omegai',[3 1]);
omegaj=sym('omegaj',[3 1]);
Rw=POCReviewRwz(wij(1),wij(2),0);
tRw=Rw(1:2,1:2);
etai=omegai(1:2);
etaj=omegaj(1:2);
omegaij=omegai-Rw*omegaj;
etaij=etai-tRw*etaj;

simplify(tRw*wij-wij)

simplify(wij.'*etaij-twij.'*omegaij)
