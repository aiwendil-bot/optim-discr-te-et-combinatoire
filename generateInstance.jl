mutable struct instance

    n ::Int64 #nb objets
    m ::Int64 #nb sacs
    c #coûts
    w #poids
    W::Vector{Int64} #capacités

end

function generateRandomlyInstanceM01KP(n = 100, m = 20, max_ci = 100, max_wi = 100, max_cap = 50)

    #verboseUtility = false # rapporte (ou pas) les items par ordre decroissant

    # --- creation de l'instance
    rnd_c = rand(1:max_ci,n); # c_i \in [1,max_ci]
    rnd_w = rand(1:max_wi,n) # w_i \in [1,max_wi]
    rnd_cap = rand(10:max_cap,m)
    sort!(rnd_cap, rev=true)
    # rank the items according the decreasing values u_i = c_i/w_i
    utilite = rnd_c ./ rnd_w
    reord = sortperm(utilite, rev=true)
    ukp   = instance(n, m, zeros(n), zeros(n), zeros(m))
    for i = 1:n
        ukp.c[i] = rnd_c[reord[i]]
        ukp.w[i] = rnd_w[reord[i]]
        #=
        if (verboseUtility == true)
            @printf "(%d %d %.2f) \n " ukp.c[i] ukp.w[i] utilite[reord[i]]
        end
        =#
    end
            
    ukp.W = rnd_cap
                
    return ukp
end
