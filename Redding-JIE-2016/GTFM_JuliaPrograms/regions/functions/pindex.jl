function pindex(param, fund, w, dtradesh, nobs)

    global alpha sigma theta epsilon LL

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H
    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]

    # gamma function
    gammaf = gamma.((theta+1-sigma)./theta)

    # price index
    P = ((gammaf.^-theta).*a.*(w.^-theta)./dtradesh).^(-1./theta)

    return P
end