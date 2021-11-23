using LinearAlgebra

function backtracking(objets)
    variables_fixées = findall(x -> x == 1, objets)
    return variables_fixées[end]
end

function branchandbound(couts::Vector{Float64},poids::Vector{Float64}, capacite::Int64)
    compteurnoeuds::Int64 = 0
    bornemax = 0
    bornemin = 0
    solrelax = zeros(length(poids))
    solglouton = zeros(length(poids))
    i = 1
    while i < length(poids)
        if poids[i] > capacite
            break
        end
        solrelax[i] = 1

        if dot(solrelax,poids) + poids[i+1] > capacite
            solrelax[i+1] = (capacite - dot(solrelax, poids)) / poids[i+1]
            bornemax = floor(dot(solrelax, couts))
            solrelax[i+1] = 0
            break
        end
        i += 1
    end
    sommepoids = 0
    k= 1
    while k <= length(poids)
        if sommepoids + poids[k] <= capacite
            compteurnoeuds += 1
            solglouton[k] = 1
            sommepoids += poids[k]
        end
        k += 1
    end
    bornemin = dot(solglouton, couts)
    objets_pris = findall(x -> x == 1, solglouton)
    soltemp = copy(solglouton)
    variable_branchement = objets_pris[end]
    explor = true

    while explor && bornemin != bornemax
        soltemp[variable_branchement] = 0
        pop!(objets_pris)
        sommepoids = dot(soltemp, poids)

        if variable_branchement == length(poids)
            variable_branchement = backtracking(soltemp)
            continue
        end
        solcalcul = copy(soltemp)
        sommepoidscalcul = copy(sommepoids)
        objetspriscalculs = copy(objets_pris)
        i = variable_branchement + 1
        while i <= length(poids)
            if sommepoidscalcul + poids[i] <= capacite
                solcalcul[i] = 1
                sommepoidscalcul += poids[i]
                push!(objetspriscalculs, i)
            else
                break
            end
            i += 1
        end
        U0 = (objetspriscalculs[end] < length(poids) - 1) ? (dot(solcalcul, couts) + floor((capacite - sommepoidscalcul)/poids[objetspriscalculs[end] + 2] * couts[objetspriscalculs[end] + 2])) : dot(solcalcul, couts)
        U1 = (objetspriscalculs[end] < length(poids)) ? (dot(solcalcul, couts) + floor(couts[objetspriscalculs[end] + 1] - (poids[objetspriscalculs[end] + 1] - (capacite - sommepoidscalcul)) * couts[objetspriscalculs[end]] / poids[objetspriscalculs[end]])) : dot(solcalcul, couts)
        borne_martello = max(U0, U1)
        println("i : ",i-1)
        println(variable_branchement)
        println(soltemp)
        println("U0")
        println(U0)
        println("U1")
        println(U1)
        println(borne_martello)
        if borne_martello <= bornemin
            if length(objets_pris) < 1
                explor = false
            else
                variable_branchement = backtracking(soltemp)
            end
            continue
        else #si c'est intéressant


            for k in (variable_branchement + 1):length(poids)  #glouton

                if sommepoids + poids[k] <= capacite
                    compteurnoeuds += 1
                    soltemp[k] = 1
                    sommepoids += poids[k]
                    push!(objets_pris, k)
                end
            end #fin glouton
            variable_branchement = objets_pris[end]
            if dot(soltemp, couts) > bornemin
                bornemin = dot(soltemp, couts)
                solglouton = copy(soltemp)
            end
        end
        if length(objets_pris) < 1
            explor = false
        end
    end

    return solglouton, dot(solglouton,couts), compteurnoeuds
end

mutable struct instance

    n ::Int64 #nb objets
    c #coûts
    w #poids
    W::Int64 #capacités

end

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
    end

    ukp.W = round(Int64, sum(ukp.w)/2)

    return ukp
end
