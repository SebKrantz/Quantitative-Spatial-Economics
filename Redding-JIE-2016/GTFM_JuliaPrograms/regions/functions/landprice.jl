function landprice(param, fund, L, w)
    # global alpha sigma theta epsilon LL

    # parameters
    alpha = param[1]
    H = fund[:,3]

    r = ((1-alpha)/alpha) .* ((w .* L) ./ H)

    return r
end