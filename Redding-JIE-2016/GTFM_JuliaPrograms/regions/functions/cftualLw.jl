
function cftualLw(param, fund, L, w, tradesh, dist, Cdist, nobs)
    # global alpha sigma theta epsilon LL

    xtic = time()

    alpha = param[1]
    theta = param[2]
    epsilon = param[3]

    a = fund[:, 1]
    b = fund[:, 2]
    H = fund[:, 3]

    dd = dist .^ (-theta)
    Cdd = Cdist .^ (-theta)
    ddhat = Cdd ./ dd

    gdp = w .* L
    dtradesh = diag(tradesh)
    lambda = L ./ LL

    Cwconverge = 0
    CLconverge = 0

    CL_i = L
    Cw_i = w

    dd = dist .^ (-theta)

    println(">>>> Start Counterfactual Wage and Population Convergence <<<<")

    xx = 1
    while xx < 2000
        x = 1
        while x < 2000
            pwmat = ((Cw_i ./ w) .^ (-theta)) * ones(1, nobs)
            nummat = tradesh .* ddhat .* pwmat
            denom = sum(nummat, dims = 1)
            denommat = ones(nobs, 1) * denom
            Ctradesh = nummat ./ denommat
            test = sum(Ctradesh, dims = 1)
            mntest = mean(test)
            income = (Cw_i ./ w) .* (CL_i ./ L) .* gdp
            expend = Ctradesh * income
            income_r = round.(income .* (10 .^ 6))
            expend_r = round.(expend .* (10 .^ 6))
            if income_r == expend_r
                x = 10000
                Cwconverge = 1
            else
                Cw_e = Cw_i .* (expend ./ income) .^ (1 / theta)
                Cw_i = 0.25 .* Cw_e + 0.75 .* Cw_i
                Cw_i = Cw_i ./ geomean(Cw_i)
                Cwconverge = 0
                x = x + 1
            end
        end
        dCtradesh = diag(Ctradesh)
        num = ((dCtradesh ./ dtradesh) .^ (-alpha * epsilon / theta)) .* ((CL_i ./ L) .^ (-epsilon * (1 - alpha))) .* lambda
        CL_e = zeros(nobs, 1)
        CL_e = num ./ sum(num)
        CL_e = CL_e .* LL
        CL_i_r = round.(CL_i .* (10 .^ 6))
        CL_e_r = round.(CL_e .* (10 .^ 6))
        if CL_i_r == CL_e_r
            xx = 10000
            CLconverge = 1
        else
            CL_e = CL_i .* (CL_e ./ CL_i) .^ (1 / (epsilon * (1 - alpha)))
            CL_i = 0.25 .* CL_e + 0.75 .* CL_i
            CLconverge = 0
            xx = xx + 1
        end
    end

    xtic = time() - xtic
    return Cw_i, CL_i, Ctradesh, dCtradesh, CLconverge, Cwconverge, xtic
end
