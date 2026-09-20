# OLS regression;
# Julia replacement for the MATLAB Statistics Toolbox regress();
# X must already contain a constant column if one is wanted;

using LinearAlgebra
using Statistics: mean
using SpecialFunctions: beta_inc

# Student-t distribution function
function tcdf(t, nu)
    x = nu / (nu + t^2)
    p = 0.5 * beta_inc(nu / 2, 0.5, x)[1]
    return t >= 0 ? 1 - p : p
end

# Student-t inverse distribution function, by bisection
function tinv(p, nu)
    lo = -1.0e3
    hi = 1.0e3
    for _ in 1:200
        mid = (lo + hi) / 2
        if tcdf(mid, nu) < p
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

# F distribution upper tail
function fccdf(f, d1, d2)
    f <= 0 && return 1.0
    x = d2 / (d2 + d1 * f)
    return beta_inc(d2 / 2, d1 / 2, x)[1]
end

# b     : coefficients
# bint  : confidence intervals for the coefficients
# r     : residuals
# rint  : intervals for the residuals (outlier diagnostic)
# stats : (R2, F, p, s2)
function regress(y, X, siglevel = 0.05)

    y = vec(y)
    n, k = size(X)

    b = X \ y
    r = y .- X * b

    dof = n - k
    s2 = dot(r, r) / dof
    XtXinv = inv(Symmetric(Matrix(X' * X)))

    tcrit = tinv(1 - siglevel / 2, dof)

    # Coefficient intervals
    se = sqrt.(max.(diag(XtXinv) .* s2, 0.0))
    bint = hcat(b .- tcrit .* se, b .+ tcrit .* se)

    # Residual intervals, using the leverage of each observation
    h = diag(X * XtXinv * X')
    rse = sqrt.(max.(s2 .* (1 .- h), 0.0))
    rint = hcat(r .- tcrit .* rse, r .+ tcrit .* rse)

    # Goodness of fit, measured about the mean when a constant is included
    hascons = any(j -> allequal(view(X, :, j)) && X[1, j] != 0, 1:k)
    sstot = hascons ? sum(abs2, y .- mean(y)) : sum(abs2, y)
    ssres = dot(r, r)
    R2 = 1 - ssres / sstot

    dfmod = hascons ? k - 1 : k
    F = dfmod > 0 ? ((sstot - ssres) / dfmod) / s2 : NaN
    p = dfmod > 0 ? fccdf(F, dfmod, dof) : NaN

    return b, bint, r, rint, (R2 = R2, F = F, p = p, s2 = s2)
end
