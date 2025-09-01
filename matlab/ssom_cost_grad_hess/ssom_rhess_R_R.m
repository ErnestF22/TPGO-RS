function hrr = ssom_rhess_R_R(R, Rdot, T, lambdas, problem_data)
    egR = ssom_egrad_R(R, T, lambdas, problem_data);
    hrr_1 = zeros(size(R));
    N = problem_data.sz(3);
    for ii = 1:N
        R_i = R(:,:,ii);
        Rdot_i = Rdot(:,:,ii);
        egRi = egR(:,:,ii);
        tmp_ii_1 = Rdot_i * R_i' * egRi + R_i * Rdot_i' * egRi;
        tmp_ii_2 = R_i * egRi' * Rdot_i + Rdot_i * egRi' * R_i;
        hrr_1(:,:,ii) = -0.5 * (tmp_ii_1 + tmp_ii_2); %!! -
    end
    hrr = hrr_1;
end