function testnetwork_utils

N = 50;
mindeg = 3;

testdata = testNetwork_params(3, N, 'banded', mindeg);

% all for base node id = 1

Rs = G2R(testdata.gi);
R1_offset = Rs(:,:,1);

Ts = G2T(testdata.gitruth);
T1_offset = Ts(:,1);

%% from gij to gi

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


%% and back (gi to gij_self, compare gij to gij_self)

% R

giR = G2R(testdata.gi);
gijR_self = from_gi_to_gij_R(giR, testdata.E);

disp("max(abs(giR - gijR_self))")
disp(max(abs(gijR - gijR_self), [], "all"))

% T

gijT_self = from_gi_to_gij_T(testdata.gi, testdata.E);

disp("max(abs(giT - gijT_self))")
disp(max(abs(gijT - gijT_self), [], "all"))

%% gi and gitruth, gij and gijtruth, lambdaij and lambdaijtruth are the same?

testdata_gijtruth_with_translations_normalized = testdata.gijtruth;

for ee = 1:size(testdata.E, 1)
    testdata_gijtruth_with_translations_normalized(1:3, end, ee) = ...
        testdata.gijtruth(1:3, end, ee) / norm(testdata.gijtruth(1:3, end, ee));
end

disp("max(abs(testdata.gij - testdata_gijtruth_with_translations_normalized), [], ""all"")")
disp(max(abs(testdata.gij - testdata_gijtruth_with_translations_normalized), [], "all"))

disp("max(abs(testdata.lambdaij - testdata.lambdaijtruth), [], ""all"")")
disp(max(abs(testdata.lambdaij - testdata.lambdaijtruth), [], "all"))

disp("max(abs(testdata.gi - testdata.gitruth), [], ""all"")")
disp(max(abs(testdata.gi - testdata.gitruth), [], "all"))

end %file function


% Notes:
%
% - testdata.gijtruth and testdata.gi are equal, but in testdata.gi the
% translations are normalized
% - to go from giT to gijT and vice-versa, rotations are needed to (the
% transformations need to undergo a similar procedure to the rotations, but
% with affine matrices, and then the translations are taken from the 
% resulting translation part)