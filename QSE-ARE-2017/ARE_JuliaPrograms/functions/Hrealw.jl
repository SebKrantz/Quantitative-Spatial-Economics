function Hrealw(fund, L, w, tradesh)
    
    # global alpha sigma tLL F

    a = fund[:, 1]
    H = fund[:, 2]

    dtradesh = diag(tradesh)

    realwage = ((L ./ (sigma .* F .* dtradesh)).^(alpha ./ (sigma - 1))) .* (a .^ alpha) .* ((L ./ H).^(-(1 - alpha)))
    realwage = realwage ./ (alpha .* ((sigma ./ (sigma - 1)).^alpha) .* (((1 - alpha) ./ alpha).^(1 - alpha)))

    return realwage
end