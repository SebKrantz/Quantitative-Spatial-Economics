function pindex(param, fund, w, dtradesh)

    # global alpha sigma theta epsilon LL

    # parameters
    theta = param[2]

    a = fund[:,1]

    # gamma function
    gammaf = gamma.((theta + 1 - sigma) / theta)

    # price index
    P = ((gammaf .^ -theta) .* a .* (w .^ -theta) ./ dtradesh) .^ (-1 / theta)

    return P
end