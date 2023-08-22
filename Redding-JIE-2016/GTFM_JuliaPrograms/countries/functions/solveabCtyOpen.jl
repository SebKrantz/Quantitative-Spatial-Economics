function solveabCtyOpen(param, observe, dwght, dist, nobs)
    global alpha sigma theta epsilon LL LLwest LLeast

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
            pwmat = a_i .* (observe[:, 2] .^ (-theta))
            nummat = dd .* pwmat
            denom = sum(nummat)
            denommat = ones(nobs) * denom
            tradesh = nummat ./ denommat

            test = sum(tradesh)
            mntest = mean(test)

            income = observe[:, 2] .* observe[:, 1]
            expend = tradesh * income

            income_r = round.(income .* (10 .^ 6))
            expend_r = round.(expend .* (10 .^ 6))

            println([x, maximum(abs.(income_r - expend_r))])

            if income_r == expend_r
                println(">>>> Productivity Convergence Achieved <<<<")
                x = 10000
                aconverge = 1
            else
                a_e = a_i .* (income ./ expend)
                a_i = 0.25 .* a_e + 0.75 .* a_i
                a_i[Iwest .== 1] = a_i[Iwest .== 1] ./ geomean(a_i[Iwest .== 1])
                aconverge = 0
                x = x + 1
            end
        end

        dtradesh = diag(tradesh)

        num = b_i .* ((a_i ./ dtradesh) .^ (alpha * epsilon / theta)) .* ((observe[:, 1] ./ observe[:, 3]) .^ (-epsilon * (1 - alpha)))
        L_e = zeros(nobs)
        L_e[Iwest .== 1] = num[Iwest .== 1] ./ sum(num[Iwest .== 1])
        L_e[Ieast .== 1] = num[Ieast .== 1] ./ sum(num[Ieast .== 1])
        L_e[Iwest .== 1] = L_e[Iwest .== 1] .* LLwest
        L_e[Ieast .== 1] = L_e[Ieast .== 1] .* LLeast

        L_r = round.(observe[:, 1] .* (10 .^ 6))
        L_e_r = round.(L_e .* (10 .^ 6))
        gap = maximum(abs.(L_e_r - L_r))

        println([xx, gap])

        if gap == 0
            println(">>>> Population Convergence Achieved <<<<")
            xx = 10000
            bconverge = 1
        else
            b_e = b_i .* (observe[:, 1] ./ L_e)
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