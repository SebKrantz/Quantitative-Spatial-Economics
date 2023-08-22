function realwage(param, fund, L, w, tradesh, dist, nobs)

    global alpha sigma theta epsilon LL

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
    # gamma function
    gammaf = gamma((theta+1-sigma)./theta)

    # domestic trade share
    dtradesh = diag(tradesh)

    # real wage
    realwage = ((a./dtradesh).^(alpha./theta)).*((L./H).^(-(1-alpha)))
    realwage = realwage./(alpha.*(gammaf.^alpha).*(((1-alpha)./alpha).^(1-alpha)))

    return realwage

end