function script_tmp

for ii = 1:1000
    [lambda, lambda_pim, eigenvalue_check_ok, imaginary_eigenvalues] = ...
        test_Hmat_ssom_proj;
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

end