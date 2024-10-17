function POCTranslationSfMPositiveScalesConstraint

u=cnormalize(randn(2,3));

cvx_begin
    variable lambda(3)
    minimize norm(u*lambda)
    subject to
        lambda>=1
cvx_end

