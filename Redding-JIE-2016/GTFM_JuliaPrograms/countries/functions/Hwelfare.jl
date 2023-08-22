function Hwelfare(param, fund, L, w, tradesh, dist, nobs)

    global alpha theta epsilon Hsigma LL LLwest LLeast F

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

    # domestic trade share
    dtradesh = diagm(tradesh)

    # welfare
    welf = deltaf.*(b.^(1./epsilon)).*(a.^alpha).*((1./dtradesh).^(alpha./theta)).*(H.^(1-alpha))
    welf = welf.*(L.^(-((1./epsilon)+(1-alpha)-(alpha./theta))))
    welf = welf./(alpha.*(((1-alpha)./alpha).^(1-alpha)).*((Hsigma./(Hsigma-1)).^alpha).*((Hsigma.*F).^(alpha./theta)))
    welf[Iwest.==1] = welf[Iwest.==1]./(LLwest.^(-1./epsilon))
    welf[Ieast.==1] = welf[Ieast.==1]./(LLeast.^(-1./epsilon))

    return welf
end