function Hlandprice(fund, L, w)
    # global alpha sigma LL LLwest LLeast

    a = fund[:, 1]
    H = fund[:, 2]

    r = ((1 - alpha) / alpha) * ((w .* L) ./ H)

    return r
end