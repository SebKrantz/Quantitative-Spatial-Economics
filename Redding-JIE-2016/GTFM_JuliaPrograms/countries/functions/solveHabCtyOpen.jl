function solveHabCtyOpen(param, observe, dwght, dist, nobs)
    global alpha theta epsilon LL LLwest LLeast

    xtic = time()

    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    aconverge = 0
    bconverge = 0

    a_i = ones(nobs)
    b_i = ones(nobs)

    dd = dist .^ (-theta)
    dd = dwght .* dd

    println(">>>> Start productivity and amenities Convergence <<<<")

    xx = 1
    while xx < 2000
        x = 1
        while x < 2000
            pwmat = L .* (a_i .^ theta) .* (w .^ (-theta)) .* ones(1, nobs)
            nummat = dd .* pwmat
            denom = sum(nummat)
            denommat = ones(nobs, 1) * denom
            tradesh = nummat ./ denommat

            test = sum(tradesh)
            mntest = mean(test)

            income = w .* L
            expend = tradesh * income

            income_r = round.(income .* (10 .^ 6))
            expend_r = round.(expend .* (10 .^ 6))

            [x, maximum(abs.(income_r - expend_r))]

            if income_r == expend_r
                x = 10000
                aconverge = 1
            else
                a_e = a_i .* ((income ./ expend) .^ (1 / theta))
                a_i = 0.25 .* a_e + 0.75 .* a_i
                a_i[Iwest .== 1] = a_i[Iwest .== 1] ./ geomean(a_i[Iwest .== 1])
                aconverge = 0
                x = x + 1
            end
        end

        dtradesh = diag(tradesh)

        num = b_i .* (a_i .^ (alpha * epsilon)) .* (H .^ (epsilon * (1 - alpha))) .* (dtradesh .^ (-alpha * epsilon / theta)) .* (L .^ (-((epsilon * (1 - alpha)) - (alpha * epsilon / theta))))
        L_e = zeros(nobs)
        L_e[Iwest .== 1] = num[Iwest .== 1] ./ sum(num[Iwest .== 1])
        L_e[Ieast .== 1] = num[Ieast .== 1] ./ sum(num[Ieast .== 1])
        L_e[Iwest .== 1] = L_e[Iwest .== 1] .* LLwest
        L_e[Ieast .== 1] = L_e[Ieast .== 1] .* LLeast

        L_r = round.(L .* (10 .^ 6))
        L_e_r = round.(L_e .* (10 .^ 6))
        gap = maximum(abs.(L_e_r - L_r))

        [xx, gap]

        if gap == 0
            xx = 10000
            bconverge = 1
        else
            b_e = b_i .* (L ./ L_e)
            b_i = 0.25 .* b_e + 0.75 .* b_i
            b_i[Iwest .== 1] = b_i[Iwest .== 1] ./ geomean(b_i[Iwest .== 1])
            b_i[Ieast .== 1] = b_i[Ieast .== 1] ./ geomean(b_i[Ieast .== 1])
            bconverge = 0
            xx = xx + 1
        end
    end

    xtic = time() - xtic
    xtic
end