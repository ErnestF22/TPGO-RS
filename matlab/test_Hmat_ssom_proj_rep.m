function test_Hmat_ssom_proj_rep

eigenvalue_valid_tests = 0;
for ii = 1:200
    [lambda, lambda_pim, v_out, v_out_pim, eigenvalue_check_ok, imaginary_eigenvalues] = ...
        test_Hmat_ssom_proj;
    if eigenvalue_check_ok
        eigenvalue_valid_tests = eigenvalue_valid_tests + 1;
        disp("abs(lambda-lambda_pim)")
        disp(abs(lambda-lambda_pim))
    else
        if lambda_pim > lambda - 1e-2
            error("lambda_pim > lambda")
        end
    end
    disp("ii in script_tmp")
    disp(ii)
    if ((~is_equal_floats(lambda_pim, lambda) && eigenvalue_check_ok) || imaginary_eigenvalues)
        disp("abs(lambda_pim - lambda)")
        disp(abs(lambda_pim - lambda))
        disp("eigenvalue_check_ok")
        disp(eigenvalue_check_ok)
        disp("imaginary_eigenvalues")
        disp(imaginary_eigenvalues)
        error("lambda_pim != lambda")
    end
end

disp("eigenvalue_valid_tests")
disp(eigenvalue_valid_tests)
end