function solveHLwCtyOpen(fund, dist, bord, bordc, nobs)

    # global alpha sigma LL LLwest LLeast

    xtic = time()

    # Extract location characteristics from fundamentals matrix
    a = fund[:, 1]
    H = fund[:, 2]
    Iwest = fund[:, 3]
    Ieast = fund[:, 4]

    # convergence indicator
    converge = 0
    tradesh = 0 
    dtradesh = 0

    # Initialization based on a symmetric allocation
    L_i = ones(nobs) .* (LL / nobs)
    w_i = ones(nobs)

    # trade costs
    dd = (dist .* bord .* bordc) .^ (1 - sigma)

    # START LOOP TO SOLVE FOR WAGES AND POPULATION
    x = 1
    while x < 2000

        # Trade share
        pwmat = (L_i .* (a .^ (sigma - 1)) .* (w_i .^ (1 - sigma))) * ones(1, nobs)
        nummat = dd .* pwmat
        denom = sum(nummat)
        denommat = ones(nobs) * denom
        tradesh = nummat ./ denommat

        # Income equals expenditure equilibrium condition
        income = w_i .* L_i
        expend = tradesh * income

        # domestic trade share
        dtradesh = diag(tradesh)

        # Population mobility equilibrium condition
        num = ((a .^ alpha) .* (H .^ (1 - alpha)) .* (dtradesh .^ (-alpha / (sigma - 1)))) .^ ((sigma - 1) / ((sigma * (1 - alpha)) - 1))
        L_e = zeros(nobs)
        L_e[Iwest .== 1] = num[Iwest .== 1] ./ sum(num[Iwest .== 1])
        L_e[Ieast .== 1] = num[Ieast .== 1] ./ sum(num[Ieast .== 1])
        L_e[Iwest .== 1] = L_e[Iwest .== 1] .* LLwest
        L_e[Ieast .== 1] = L_e[Ieast .== 1] .* LLeast

        # Convergence criterion
        income_r = round.(income .* (10 .^ 6))
        expend_r = round.(expend .* (10 .^ 6))
        L_i_r = round.(L_i .* (10 .^ 6))
        L_e_r = round.(L_e .* (10 .^ 6))

        [x, maximum(abs.(income_r - expend_r)), maximum(abs.(L_e_r - L_i_r))]

        # Update loop
        if income_r == expend_r && L_i_r == L_e_r
            # display('>>>> Convergence Achieved <<<<')
            x = 10000
            converge = 1
        else
            w_e = w_i .* (expend ./ income) .^ (1 / (sigma - 1))
            w_i = (0.25 .* w_e) + (0.75 .* w_i)
            L_i = (0.25 .* L_e) + (0.75 .* L_i)
            # Normalization
            # Choose geometric mean wage in West as numeraire
            w_i[Iwest .== 1] = w_i[Iwest .== 1] ./ geomean(w_i[Iwest .== 1])
            wconverge = 0
            x = x + 1
        end

    end

    # END LOOP TO SOLVE FOR POPULATION AND WAGES

    xtic = time() - xtic

    return w_i, L_i, tradesh, dtradesh, converge, xtic

end