using LinearAlgebra

function backtracking(objets)
    variables_fixees = findall(x -> x == 1, objets)
    return variables_fixees[end]

end

function solver_O1UKP_V1(couts::Vector{Float64},poids::Vector{Float64}, capacite::Int64)
    compteurnoeuds::Int64 = 0
    bornemin::Int64 = 0
    solrelax = zeros(length(poids))
    solglouton = zeros(length(poids))
    s::Int64 = calcul_dernier_objet_non_bloquant(couts, poids, capacite)
    c_barre = capacite - sum([poids[i] for i in 1:(s-1)])
    bornemax::Int64 = sum([couts[i] for i in 1:(s-1)]) +  floor(c_barre * couts[s]/poids[s])
    println("borne max : ", bornemax)
    sommepoids = 0
    k= 1
    while k <= length(poids)
        if sommepoids + poids[k] <= capacite
            compteurnoeuds +=1
            solglouton[k] = 1
            sommepoids += poids[k]
        end
        k += 1
    end
    bornemin = dot(solglouton, couts)
    objets_pris = findall(x -> x == 1, solglouton)
    soltemp = copy(solglouton)
    indice_variables = [i for i in 1:length(poids)]
    variable_branchement = objets_pris[end]
    explor = true
    while explor && bornemin != bornemax
        soltemp[variable_branchement] = 0
        if length(objets_pris) >0
            pop!(objets_pris)
        end
        s = variable_branchement
        l = 1
        sommepoids = variable_branchement ==1 ? 0 : sum([poids[i]*soltemp[i] for i in 1:(variable_branchement - 1)])
        for k in (variable_branchement + 1):length(soltemp)
            if sommepoids + poids[k] <= capacite
                s = k
                sommepoids += poids[k]
                k += 1
            else
                break
            end
        end
        borne_LPK = 0

        if s == length(soltemp)
            borne_LPK = dot(soltemp, couts) + couts[s]
        else
            capacite_residuelle = capacite - sum([poids[i] for i in findall(x->x==1, soltemp)])
            z_current = sum([couts[i] for i in findall(x->x==1,soltemp)])
            alpha = sum([couts[j] for j in (variable_branchement + 1):s]) + (capacite_residuelle - sum([poids[j] for j in (variable_branchement + 1):s])) / poids[s + 1] * couts[s + 1]
            borne_LPK = z_current + alpha
        end

        if borne_LPK <= bornemin
            if length(objets_pris) < 1
                explor = false
            else
                variable_branchement = backtracking(soltemp)
            end
            continue
        else
            sommepoidsb = variable_branchement == 1 ? 0 : sum([poids[i]*soltemp[i] for i in 1:(variable_branchement - 1)])
            for k in (variable_branchement + 1):length(poids)  #glouton

                if sommepoidsb + poids[k] <= capacite
                    compteurnoeuds += 1
                    soltemp[k] = 1
                    sommepoidsb += poids[k]
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
    return dot(solglouton,couts), compteurnoeuds
end

function calcul_dernier_objet_non_bloquant(couts::Vector{Float64},poids::Vector{Float64}, capacite::Int64)::Int64
    s = 1
    sommepoids::Int64 = 0
    while sommepoids + poids[s] <= capacite && s < length(poids)
        sommepoids += poids[s]
        s += 1
    end
    return s
end

mutable struct instance

    n ::Int64 #nb objets
    c #coûts
    w #poids
    W :: Int64

end

function generateRandomlyInstance01UKP(n = 100, max_ci = 100, max_wi = 100)

    #verboseUtility = false # rapporte (ou pas) les items par ordre decroissant

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
        #=
        if (verboseUtility == true)
            @printf "(%d %d %.2f) \n " ukp.c[i] ukp.w[i] utilite[reord[i]]
        end
        =#
    end

    ukp.W = round(Int64, sum(ukp.w)/2)

    return ukp
end
