function Hwelfaregains(ctradesh, tradesh, cL, L)

    # global alpha sigma LL

    # domestic trade share
    dtradesh = diag(tradesh)
    cdtradesh = diag(ctradesh)

    # welfare gains
    welfgain = ((dtradesh ./ cdtradesh) .^ (alpha / (sigma - 1))) .* ((L ./ cL) .^ (((sigma * (1 - alpha)) - 1) / (sigma - 1)))

    return welfgain
end