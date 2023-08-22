function solveImmobileCtyOpen(param, fund, L, dwght, dist, nobs)

    global alpha sigma theta epsilon LL

    xtic = time()

    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    a = fund[:, 1]
    b = fund[:, 2]
    H = fund[:, 3]
    Iwest = fund[:, 4]
    Ieast = fund[:, 5]

    wconverge = 0

    w_i = ones(nobs)

    dd = dist .^ (-theta)
    dd = dwght .* dd

    println(">>>> Start Wage and Population Convergence <<<<")

    x = 1
    while x < 2000

        pwmat = a .* (w_i .^ (-theta)) * ones(1, nobs)
        nummat = dd .* pwmat
        denom = sum(nummat)
        denommat = ones(nobs) * denom
        tradesh = nummat ./ denommat

        test = sum(tradesh)
        mntest = mean(test)

        income = w_i .* L
        expend = tradesh * income

        income_r = round.(income .* (10 .^ 6))
        expend_r = round.(expend .* (10 .^ 6))

        println([x, maximum(abs.(income_r - expend_r))])

        if income_r == expend_r
            x = 10000
            wconverge = 1
        else
            w_e = w_i .* (expend ./ income) .^ (1 / theta)
            w_i = 0.25 .* w_e + 0.75 .* w_i
            w_i = w_i ./ geomean(w_i[Iwest .== 1])
            wconverge = 0
            x = x + 1
        end

    end

    dtradesh = diag(tradesh)

    gammaf = gamma((theta + 1 - sigma) / theta)

    P = ((gammaf .^ -theta) .* a .* (w_i .^ -theta) ./ dtradesh) .^ (-1 / theta)

    r = ((1 - alpha) ./ alpha) .* ((w_i .* L) ./ H)

    xtic = time() - xtic
    return w_i, P, r, tradesh, dtradesh, wconverge, xtic
end