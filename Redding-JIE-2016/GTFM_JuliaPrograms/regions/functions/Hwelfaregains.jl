function Hwelfaregains(param, Ctradesh, tradesh, CL, L)

    # global alpha Hsigma theta epsilon LL

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # domestic trade share
    dtradesh = diag(tradesh)
    Cdtradesh = diag(Ctradesh)

    # welfare gains
    welfgain = ((dtradesh ./ Cdtradesh) .^ (alpha/theta)) .* ((L ./ CL) .^ ((1/epsilon)+(1-alpha)-(alpha/theta)))
    return welfgain
end