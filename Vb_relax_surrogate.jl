using LinearAlgebra

function backtracking(objets, F)
    variables_fixeesb = findall(x -> x == 1, objets)
    variables_fixees = intersect(variables_fixeesb, F)
    return length(variables_fixees) > 0 ? variables_fixees[end] : 0
end

function branchandbound(couts::Vector{Float64},poids::Vector{Float64}, capacite::Float64, F::Vector{Int64})
    compteurnoeuds::Int64 = 0
    bornemax = 0
    bornemin = 0
    solrelax = ones(length(poids))
    solglouton = ones(length(poids))
    for idx in F
        solrelax[idx], solglouton[idx] = 0,0
    end
    s = 1
    alpha::Float64 = 0
    sommepoids = dot(solrelax, poids)
    #=
    s+=1
    c_barre = capacite - sum([poids[i] for i in 1:(s-1)])
    U_0 = sum([couts[i] for i in 1:(s-1)]) + floor(c_barre * couts[s + 1]/poids[s + 1])
    U_1 = sum([couts[i] for i in 1:(s-1)]) + floor(couts[s] - (poids[s] - c_barre) * couts[s - 1] / poids[s - 1] )
    bornemax = max(U_0, U_1)
    =#

    for k in F
        if sommepoids + poids[k] <= capacite
            compteurnoeuds +=1
            solglouton[k] = 1
            sommepoids += poids[k]
        end
    end
    bornemin = dot(solglouton, couts)
    objets_pris = findall(x -> x == 1, solglouton)
    soltemp = copy(solglouton)
    variable_branchement = intersect(objets_pris, F)[end]
    explor = true
    while explor
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
                variable_branchement = backtracking(soltemp, F)
                if variable_branchement == 0
                    explor = false
                end
            end
            continue
        else
            sommepoidsb = dot(poids,soltemp)
            for k in intersect(F,(variable_branchement + 1):length(poids))  #glouton
                if sommepoidsb + poids[k] <= capacite
                    compteurnoeuds += 1
                    soltemp[k] = 1
                    sommepoidsb += poids[k]
                    deleteat!(F, findall(x->x==k, F))
                    push!(objets_pris, k)
                end
            end #fin glouton
            sort!(F)
            if length(intersect(F,objets_pris)) > 0
                variable_branchement = intersect(F,objets_pris)[end]
            else
                break
            end
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


function borneduale_surrogate(couts::Vector{Float64}, poids::Vector{Float64}, capa::Vector{Int64}, coeff::Float64, F::Vector{Int64}=Int64[])
    return branchandbound(couts, coeff .* poids, sum(coeff .* capa),F)[2]
end

