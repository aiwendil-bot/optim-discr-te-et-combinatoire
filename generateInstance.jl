function generateRandomlyInstanceUKP(n = 100, max_ci = 100, max_wi = 30)

    verboseUtility = false # rapporte (ou pas) les items par ordre decroissant

    # --- creation de l'instance
    rnd_c = rand(1:max_ci,n); # c_i \in [1,max_ci]
    rnd_w = rand(1:max_wi,n) # w_i \in [1,max_wi]
    
    # rank the items according the decreasing values u_i = c_i/w_i
    utilite = rnd_c ./ rnd_w
    reord = sortperm(utilite, rev=true)
    ukp   = instance(n, zeros(n), zeros(n), 0)
    for i = 1:n
        ukp.c[i] = rnd_c[reord[i]]
        ukp.w[i] = rnd_w[reord[i]]
        if (verboseUtility == true)
            @printf "(%d %d %.2f) \n " ukp.c[i] ukp.w[i] utilite[reord[i]]
        end
    end
            
    ukp.W = round(Int64, sum(ukp.w)/2)
                
    return ukp
end