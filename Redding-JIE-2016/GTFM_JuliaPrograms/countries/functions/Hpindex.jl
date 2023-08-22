function Hpindex(param, fund, L, w, dtradesh, nobs)
    global alpha Hsigma theta epsilon LL F

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H;
    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]

    # price index
    P = (Hsigma / (Hsigma - 1)) .* (w ./ a) .* ((L ./ (Hsigma .* F .* dtradesh)).^(-1 ./ theta))
end