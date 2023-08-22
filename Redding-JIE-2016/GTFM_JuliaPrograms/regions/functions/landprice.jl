function landprice(param, fund, L, w, d, dist, nobs)
    global alpha sigma theta epsilon LL

    xtic = tic()

    # parameters
    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    # fund[:,1]=a; fund[:,2]=b; fund[:,3]=H
    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]

    r = ((1-alpha)/alpha) * ((w*L) / H)

    return r
end