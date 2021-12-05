using LinearAlgebra

function backtracking(objets)
    variables_fixees = findall(x -> x == 1, objets)
    return variables_fixees[end]
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
    j= 1
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
        borne_dantzig = objetspriscalculs[end] < length(poids) ? (dot(solcalcul, couts) + (capacite - sommepoidscalcul)/poids[objetspriscalculs[end] + 1] * couts[objetspriscalculs[end] + 1]) : dot(solcalcul, couts)
        if borne_dantzig <= bornemin
            if length(objets_pris) < 1
                explor = false
            else
                variable_branchement = backtracking(soltemp)
            end
            continue
        else

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
branchandbound([4.0,9.0,10.0,9.0,3.0,2.0],[1.0,3.0,4.0,4.0,2.0,3.0],7)
