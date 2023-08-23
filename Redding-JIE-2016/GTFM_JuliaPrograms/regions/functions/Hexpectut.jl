function Hexpectut(param, fund, w, P, r)
    # global alpha theta epsilon LL

    # parameters
    alpha = param[1]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H
    b = fund[:,2]

    # delta function
    deltaf = gamma((epsilon-1)/epsilon)

    # expected utility
    EU = b .* (P .^ (-alpha * epsilon)) .* (r .^ (-(1-alpha) * epsilon)) .* ((w ./ alpha) .^ epsilon)
    EU = deltaf .* (sum(EU) ^ (1 / epsilon))

    return EU
end