function Hpindex(fund, L, w, dtradesh)
    
    # global alpha sigma LL F

    # fund[:,1]=a; fund[:,2]=H;
    a = fund[:,1]
    H = fund[:,2]

    # price index
    P = (sigma / (sigma - 1)) .* (w ./ a) .* ((L ./ (sigma .* F .* dtradesh)) .^ (1 / (1 - sigma)))

    return P
end