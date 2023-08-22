function solveLwCtyOpen(param,fund,dwght,dist,nobs)

    global alpha sigma theta epsilon LL LLwest LLeast

    xtic = time()

    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    a = fund[:,1]
    b = fund[:,2]
    H = fund[:,3]
    Iwest = fund[:,4]
    Ieast = fund[:,5]

    wconverge = 0
    Lconverge = 0

    L_i = ones(nobs) .* (LL/nobs)
    w_i = ones(nobs)

    dd = dist .^ (-theta)
    dd = dwght .* dd

    println(">>>> Start Wage and Population Convergence <<<<")

    xx = 1
    while xx < 2000

        x = 1
        while x < 2000

            pwmat = a .* (w_i .^ (-theta)) * ones(1, nobs)
            nummat = dd .* pwmat
            denom = sum(nummat)
            denommat = ones(nobs, 1) * denom
            tradesh = nummat ./ denommat

            test = sum(tradesh)
            mntest = mean(test)

            income = w_i .* L_i
            expend = tradesh * income

            income_r = round.(income .* (10 .^ 6))
            expend_r = round.(expend .* (10 .^ 6))

            println([x, maximum(abs.(income_r - expend_r))])

            if income_r == expend_r
                println(">>>> Wage Convergence Achieved <<<<")
                x = 10000
                wconverge = 1
            else
                w_e = w_i .* (expend ./ income) .^ (1 / theta)
                w_i = 0.25 .* w_e + 0.75 .* w_i
                w_i[Iwest .== 1] = w_i[Iwest .== 1] ./ geomean(w_i[Iwest .== 1])
                wconverge = 0
                x = x + 1
            end
        end

        dtradesh = diag(tradesh)

        num = b .* ((a ./ dtradesh) .^ (alpha * epsilon / theta)) .* ((L_i ./ H) .^ (-epsilon * (1 - alpha)))
        L_e = zeros(nobs)
        L_e[Iwest .== 1] = num[Iwest .== 1] ./ sum(num[Iwest .== 1])
        L_e[Ieast .== 1] = num[Ieast .== 1] ./ sum(num[Ieast .== 1])
        L_e[Iwest .== 1] = L_e[Iwest .== 1] .* LLwest
        L_e[Ieast .== 1] = L_e[Ieast .== 1] .* LLeast

        L_i_r = round.(L_i .* (10 .^ 6))
        L_e_r = round.(L_e .* (10 .^ 6))

        println([xx, maximum(abs.(L_e_r - L_i_r))])

        if L_i_r == L_e_r
            println(">>>> Population Convergence Achieved <<<<")
            xx = 10000
            Lconverge = 1
        else
            L_e = L_i .* (L_e ./ L_i) .^ (1 / (epsilon * (1 - alpha)))
            L_i = 0.25 .* L_e + 0.75 .* L_i
            Lconverge = 0
            xx = xx + 1
        end
    end

    xtic = time() - xtic
    println(xtic)

end