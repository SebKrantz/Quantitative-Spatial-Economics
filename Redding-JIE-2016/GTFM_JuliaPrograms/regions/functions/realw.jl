
function realw(param, fund, L, tradesh)
    # global alpha, sigma, theta, epsilon, LL

    # parameters
    alpha = param[1]
    theta = param[2]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H
    a = fund[:,1]
    H = fund[:,3]

    # gamma function
    gammaf = gamma((theta+1-sigma)/theta)

    # domestic trade share
    dtradesh = diag(tradesh)

    # real wage
    realwage = ((a ./ dtradesh) .^ (alpha/theta)) .* ((L ./ H) .^ (-(1-alpha)))
    realwage = realwage ./ (alpha .* (gammaf .^ alpha) .* (((1-alpha)/alpha)^(1-alpha)))

    return realwage
end
