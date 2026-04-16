function testnetwork_utils

N = 50;
mindeg = 3;

testdata = testNetwork_params(3, N, 'banded', mindeg); 

% for base node id = 1

Rs = G2R(testdata.gi);
R1_offset = Rs(:,:,1);

Ts = G2T(testdata.gitruth);
T1_offset = Ts(:,1);

% from gij to gi

% R
gijR = G2R(testdata.gij);
giR_self = from_gij_to_gi_R(gijR, testdata.E, testdata.NNodes, R1_offset);

% disp("[giR, giR_self]")
% disp([Rs, giR_self])

disp("max(abs(Rs - giR_self))")
disp(max(abs(Rs - giR_self), [], "all"))


% T
gijT = G2T(testdata.gijtruth);
giT_self = from_gij_to_gi_T(gijT, gijR, testdata.E, testdata.NNodes, T1_offset, R1_offset);

% disp("[giT; giT_self]")
% disp([Ts; giT_self])

disp("max(abs(Ts - giT_self))")
disp(max(abs(Ts - giT_self), [], "all"))


% and back (gi to gij_self, compare gij to gij_self)

end