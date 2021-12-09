#=
Vb_relax_surrogate:
- Julia version: 1.6.0
- Author: adrien
- Date: 2021-12-09
=#

using LinearAlgebra

function backtracking(objets, S)
    variables_fixees = findall(x -> x == 1, objets)
    for k in 1:length(variables_fixees)
        if k > length(variables_fixees)
            break
        end
        if variables_fixees[k] in S
             deleteat!(variables_fixees, k)
        end
    end
    return length(variables_fixees) > 0 ? variables_fixees[end] : 0
end

function branchandbound(couts::Vector{Float64},poids::Vector{Float64}, capacite::Float64, S::Vector{Int64})
    compteurnoeuds::Int64 = 0;
    bornemax = 0
    bornemin = 0
    solrelax = zeros(length(poids))
    solglouton = zeros(length(poids))
    for idx in S
        solrelax[idx], solglouton[idx] = 1,1
    end
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

    while explor && bornemin != bornemax && variable_branchement != 0
        soltemp[variable_branchement] = 0
        pop!(objets_pris)
        sommepoids = dot(soltemp, poids)

        if variable_branchement == length(poids)
            variable_branchement = backtracking(soltemp, S)
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

        if borne_martello <= bornemin
            if length(objets_pris) < 1
                explor = false
            else
                variable_branchement = backtracking(soltemp, S)
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



function borneduale_surrogate(couts::Vector{Float64}, poids::Vector{Float64}, capa::Vector{Int64}, coeff::Float64, S::Vector{Int64}=Int64[])
    return branchandbound(couts, coeff .* poids, sum(coeff .* capa),S)[2]
end
