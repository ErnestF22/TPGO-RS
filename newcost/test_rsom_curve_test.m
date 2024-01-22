function test_rsom_curve_test()
disp('Check that tangent is always in tangent space')
problem=test_problem();
[curve, Y0, R1, dR1, ddR1, R2, dR2, ddR2] = test_problem_curve(problem);

% Y0 = eye3d(4,3,5);

figure(1)
f=@(t) sum(multitrace(multiprod(multitransp(curve.c(t)),curve.dc(t))));
funPlot(f)
title('Should always be near zero')

figure(2);
disp('Check the first derivative: R1')
R1_st = @(t) matStack(multiprod(R1(t), Y0));
dR1_st = @(t) matStack(multiprod(dR1(t), Y0));
% funCheckDer(R1_st, dR1_st)
figure(3);
disp('Check the first derivative: R2')
R2_st = @(t) matStack(multiprod(Y0, R2(t)));
dR2_st = @(t) matStack(multiprod(Y0, dR2(t)));
% funCheckDer(R2_st, dR2_st)


figure(4);
disp('Check the second derivative: R1')
ddR1_st = @(t) matStack(multiprod(ddR1(t), Y0));
% funCheckDer(dR1_st, ddR1_st)

figure(5);
disp('Check the second derivative: R2')
ddR2_st = @(t) matStack(multiprod(Y0, ddR2(t)));
% funCheckDer(dR2_st, ddR2_st)


figure(6)
disp('Check the first derivative: full cost function')
R_st = @(t) matStack(multiprod3(R1(t), Y0, R2(t)));
dR_st = @(t) matStack(multiprod3(dR1(t), Y0, R2(t)) + multiprod3(R1(t), Y0, dR2(t)));
% funCheckDer(R_st, dR_st)

figure(7)
disp('Check the second derivative: full cost function')
ddR_st = @(t) matStack( multiprod3(ddR1(t), Y0, R2(t)) + ...
    2.*multiprod3(dR1(t), Y0, dR2(t)) + ...
    multiprod3(R1(t), Y0, ddR2(t)) );
% funCheckDer(dR_st, ddR_st)


figure(8)
disp('Check that the second derivative (Euclidean acceleration) is consistent with tangent')
dc=dcSubset(curve,[1;2]);
ddc=ddcSubset(curve,[1;2]);
funCheckDer(dc, ddc, 'angle')



figure(9)
curve_dc_t = @(t) matStack(curve.dc(t));
curve_ddc_t = @(t) matStack(curve.ddc(t));
funCheckDer(curve_dc_t, curve_ddc_t)

end %file function

function h=dcSubset(curve,idx)
    %Return handle to a function similar to curve.dc, but limited to the
    %specified indeces
    function v=dc(t)
        v=curve.dc(t);
        v=v(idx);
    end
    h=@dc;
end

function h=ddcSubset(curve,idx)
    %Return handle to a function similar to curve.ddc, but limited to the
    %specified indeces
    function v=ddc(t)
        v=curve.ddc(t);
        v=v(idx);
    end
    h=@ddc;
end