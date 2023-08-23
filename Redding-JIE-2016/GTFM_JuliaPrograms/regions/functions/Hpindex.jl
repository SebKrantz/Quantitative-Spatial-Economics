function Hpindex(param, fund, L, w, dtradesh)
    # global alpha Hsigma theta epsilon LL F

    # parameters
    theta = param[2]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H;
    a = fund[:,1]

    # price index
    P = (Hsigma / (Hsigma - 1)) .* (w ./ a) .* ((L ./ (Hsigma .* F .* dtradesh)) .^ (-1 / theta))
    return P
end