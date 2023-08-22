function landprice(param, fund, L, w, d, dist, nobs)
    global alpha sigma theta epsilon LL LLwest LLeast

    xtic = tic()

    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    a = fund[:, 1]
    b = fund[:, 2]
    H = fund[:, 3]

    r = ((1 - alpha) / alpha) * ((w * L) / H)

    return r
end