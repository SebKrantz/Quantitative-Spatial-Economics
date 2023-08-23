function expectut(param, fund, L, tradesh)

    # global alpha sigma theta epsilon LL

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H
    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]

    # delta function
    deltaf = gamma((epsilon-1)/epsilon)
    # gamma function
    gammaf = gamma((theta+1-sigma)/theta)

    # domestic trade share
    dtradesh = diag(tradesh)

    # expected utility
    EU = b .* (gammaf.^(-alpha*epsilon)) .* (alpha.^(-epsilon)) .* (((1-alpha)./alpha).^(-epsilon*(1-alpha)))
    EU = ((a./dtradesh).^(alpha*epsilon/theta)) .* ((L./H).^(-epsilon*(1-alpha))) .* EU
    EU = deltaf .* (sum(EU).^(1/epsilon))

    return EU

end