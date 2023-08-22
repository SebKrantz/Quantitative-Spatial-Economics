function welfare(param, fund, L, w, tradesh, dist, nobs)

    global alpha sigma theta epsilon LL

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H; 
    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]

    # delta function
    deltaf = gamma((epsilon-1)./epsilon)
    # gamma function
    gammaf = gamma((theta+1-sigma)./theta)

    # domestic trade share
    dtradesh = diag(tradesh)

    # welfare
    welf = deltaf.*(b.^(1./epsilon)).*((a./dtradesh).^(alpha./theta)).*(H.^(1-alpha)).*(L.^(-((1./epsilon)+(1-alpha))))
    welf = welf./(alpha.*(((1-alpha)./alpha)^(1-alpha)).*(gammaf.^alpha).*(LL.^(-1./epsilon)))

    return welf
end