function Hrealw(param, fund, L, tradesh)

    # global alpha Hsigma theta epsilon LL F

    # parameters
    alpha = param[1]
    theta = param[2]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H
    a = fund[:,1]
    H = fund[:,3]

    # domestic trade share
    dtradesh = diag(tradesh)

    # real wage
    realwage = ((L ./ (Hsigma .* F .* dtradesh)) .^ (alpha/theta)) .* (a .^ alpha) .* ((L ./ H) .^ (-(1-alpha)))
    realwage = realwage ./ (alpha*(((1-alpha)/alpha)^(1-alpha)))

    return realwage
end