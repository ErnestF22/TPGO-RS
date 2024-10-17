function rotLocFrobCostPair_gradMatrix_test
Ri=rot_randn();
Rj=rot_randn();
Rij=rot_randn();

gradc=rotLocFrobCostPair_grad(Ri,Rj,Rij);

gradMat=rotLocFrobCostPair_gradMatrix(Ri,Rj,Rij);

gradcB=multitransp(matUnstack(gradMat*[Ri';Rj']));

disp([gradc-gradcB])
