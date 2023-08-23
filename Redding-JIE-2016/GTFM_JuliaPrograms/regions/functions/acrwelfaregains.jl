function acrwelfaregains(param, Ctradesh, tradesh)
    # global alpha sigma theta epsilon LL

    # parameters
    alpha = param[1]
    theta = param[2]

    # domestic trade share
    dtradesh = diag(tradesh)
    Cdtradesh = diag(Ctradesh)

    # welfare gains
    acrwelfgain = (dtradesh ./ Cdtradesh) .^ (alpha / theta)

    return acrwelfgain
end