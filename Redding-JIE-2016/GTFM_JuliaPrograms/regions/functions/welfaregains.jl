function welfaregains(param, fund, Ctradesh, tradesh, CL, L, nobs)

    global alpha sigma theta epsilon LL

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
    # gamma function
    gammaf = gamma((theta+1-sigma)./theta)

    # domestic trade share
    dtradesh = diag(tradesh)
    Cdtradesh = diag(Ctradesh)

    # welfare gains
    welfgain = ((dtradesh./Cdtradesh).^(alpha./theta)).*((L./CL).^((1./epsilon)+(1-alpha)))

    return welfgain
end