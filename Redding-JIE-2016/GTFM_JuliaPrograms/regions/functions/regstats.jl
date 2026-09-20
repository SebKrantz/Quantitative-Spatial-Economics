# Regression diagnostics;
# Julia replacement for the MATLAB Statistics Toolbox regstats();
# x holds the regressors WITHOUT a constant; the "linear" model prepends one;
# Returns beta and covb, so that sqrt.(diag(stats.covb)) are the standard errors;

using LinearAlgebra

function regstats(y, x, model = "linear", whichstats = ["beta", "covb"])

    model == "linear" || error("regstats: only the \"linear\" model is implemented, got \"$model\"")

    y = vec(y)
    Xr = x isa AbstractVector ? reshape(x, :, 1) : x
    X = hcat(ones(eltype(float.(Xr)), size(Xr, 1)), Xr)

    n, k = size(X)

    beta = X \ y
    r = y .- X * beta
    s2 = dot(r, r) / (n - k)
    covb = s2 .* inv(Symmetric(Matrix(X' * X)))

    return (beta = beta, covb = Matrix(covb), r = r, mse = s2)
end
