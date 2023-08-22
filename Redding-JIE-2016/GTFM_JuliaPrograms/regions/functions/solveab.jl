function solveab(param, observe, dist, nobs)
    global alpha sigma theta epsilon LL

    xtic = time()

    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    L = observe[:, 1]
    w = observe[:, 2]
    H = observe[:, 3]

    aconverge = 0
    bconverge = 0

    a_i = ones(nobs)
    b_i = ones(nobs)

    dd = dist .^ (-theta)

    println(">>>> Start productivity and amenities Convergence <<<<")

    xx = 1
    while xx < 2000
        x = 1
        while x < 2000
            pwmat = a_i .* (w .^ (-theta)) .* ones(1, nobs)
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

            println([x, maximum(abs.(income_r - expend_r))])

            if income_r == expend_r
                println(">>>> Productivity Convergence Achieved <<<<")
                x = 10000
                aconverge = 1
            else
                a_e = a_i .* (income ./ expend)
                a_i = 0.25 .* a_e + 0.75 .* a_i
                a_i = a_i ./ geomean(a_i)
                aconverge = 0
                x = x + 1
            end
        end

        dtradesh = diag(tradesh)

        num = b_i .* ((a_i ./ dtradesh) .^ (alpha * epsilon / theta)) .* ((L ./ H) .^ (-epsilon * (1 - alpha)))
        L_e = zeros(nobs)
        L_e = num ./ sum(num)
        L_e = L_e .* LL

        L_r = round.(L .* (10 .^ 6))
        L_e_r = round.(L_e .* (10 .^ 6))

        println([xx, maximum(abs.(L_e_r - L_r))])

        if L_r == L_e_r
            println(">>>> Population Convergence Achieved <<<<")
            xx = 10000
            bconverge = 1
        else
            b_e = b_i .* (L ./ L_e)
            b_i = 0.25 .* b_e + 0.75 .* b_i
            b_i = b_i ./ geomean(b_i)
            bconverge = 0
            xx = xx + 1
        end
    end

    xtic = time() - xtic
    return a_i, b_i, tradesh, aconverge, bconverge, xtic
end