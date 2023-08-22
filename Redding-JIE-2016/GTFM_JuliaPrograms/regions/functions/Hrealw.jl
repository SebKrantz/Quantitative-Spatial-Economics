function Hrealw(param, fund, L, w, tradesh, dist, nobs)

    global alpha Hsigma theta epsilon LL F

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

    # domestic trade share
    dtradesh = diag(tradesh)

    # real wage
    realwage = ((L./(Hsigma.*F.*dtradesh)).^(alpha./theta)).*(a.^alpha).*((L./H).^(-(1-alpha)))
    realwage = realwage./(alpha.*(((1-alpha)./alpha).^(1-alpha)))

    return realwage

end