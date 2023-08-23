
function mobwelfaregains(param, Ctradesh, tradesh, CL, L)

    # global alpha, sigma, theta, epsilon, LL

    # parameters
    alpha = param[1]  
    theta = param[2]

    # domestic trade share 
    dtradesh = diag(tradesh)
    Cdtradesh = diag(Ctradesh)

    # welfare gains
    welfgain = ((dtradesh ./ Cdtradesh).^(alpha ./ theta)) .* ((L ./ CL).^(1 - alpha))

    return welfgain
end