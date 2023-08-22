function thetaepsopt(BETA)

    global alpha sigma theta epsilon LL nobs W
    global observe Cobserve dist0 dist1

    theta = BETA[1]
    epsilon = BETA[2]

    estparam = [alpha theta epsilon]

    # Solve for baseline region productivities and amenities
    a_i, b_i, abtradesh, aconverge, bconverge, xtic = solveab(estparam, observe, dist0, nobs)
    println(">>>> Productivity and Amenity System Converged <<<<")
    println(">>>> Check Productivity and Amenity Convergence <<<<")
    println(aconverge, " ", bconverge)
    println(">>>> Elapsed Time in Seconds <<<<")
    println(xtic)

    # Solve for counterfactual region productivities and amenities
    Ca_i, Cb_i, Cabtradesh, Caconverge, Cbconverge, Cxtic = solveab(estparam, Cobserve, dist1, nobs)
    println(">>>> Productivity and Amenity System Converged <<<<")
    println(">>>> Check Productivity and Amenity Convergence <<<<")
    println(aconverge, " ", bconverge)
    println(">>>> Elapsed Time in Seconds <<<<")
    println(xtic)

    # Change in fundamentals
    da_i = log.(Ca_i ./ a_i)
    db_i = log.(Cb_i ./ b_i)

    # Objective function
    ft = zeros(2)
    ft[1] = sum(da_i .^ 2)
    ft[2] = sum(db_i .^ 2)

    f = ft' * W * ft

    return f
end