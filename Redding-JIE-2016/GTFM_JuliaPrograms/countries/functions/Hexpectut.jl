function Hexpectut(param, fund, L, w, P, r, dist, nobs)

    global alpha theta epsilon LL

    xtic = tic()

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H; fund[:,4]=Iwest; fund[:,5]=Ieast;
    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]
    Iwest = fund[:,4]
    Ieast = fund[:,5]

    # delta function
    deltaf = gamma((epsilon-1)./epsilon)

    # expected utility
    EU = b.*(P.^(-alpha.*epsilon)).*(r.^(-(1-alpha).*epsilon)).*((w./alpha).^epsilon)
    EU[Iwest .== 1] = deltaf.*(sum(EU[Iwest .== 1]).^(1./epsilon))
    EU[Ieast .== 1] = deltaf.*(sum(EU[Ieast .== 1]).^(1./epsilon))

end