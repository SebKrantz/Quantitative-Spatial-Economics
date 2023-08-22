function Hexpectut(param, fund, L, w, P, r, dist, nobs)

    global alpha theta epsilon LL

    xtic = tic()

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H
    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]

    # delta function
    deltaf = gamma((epsilon-1)./epsilon)

    # expected utility
    EU = b.*(P.^(-alpha.*epsilon)).*(r.^(-(1-alpha).*epsilon)).*((w./alpha).^epsilon)
    EU = deltaf.*(sum(EU).^(1./epsilon))

    return EU

end