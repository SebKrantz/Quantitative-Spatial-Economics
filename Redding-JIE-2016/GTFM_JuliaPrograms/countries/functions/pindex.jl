function pindex(param, fund, w, dtradesh, nobs)
    global alpha sigma theta epsilon LL LLwest LLeast

    xtic = tic()

    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    a = fund[:, 1]
    b = fund[:, 2]
    H = fund[:, 3]

    gammaf = gamma((theta + 1 - sigma) / theta)

    P = ((gammaf .^ -theta) .* a .* (w .^ -theta) ./ dtradesh) .^ (-1 / theta)

    return P
end