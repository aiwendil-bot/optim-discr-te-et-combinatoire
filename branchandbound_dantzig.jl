using LinearAlgebra

function backtracking(objets)
    variables_fixées = findall(x -> x == 1, objets)
    return variables_fixées[end]
end

function branchandbound(couts::Vector{Int64},poids::Vector{Int64}, capacite::Int64)
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
        println(objetspriscalculs)
        println(variable_branchement)
        println(soltemp)
        borne_dantzig = objetspriscalculs[end] < length(poids) ? (dot(solcalcul, couts) + (capacite - sommepoidscalcul)/poids[objetspriscalculs[end] + 1] * couts[objetspriscalculs[end] + 1]) : dot(solcalcul, couts)
        println(borne_dantzig)
        println(objets_pris)
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
        j += 1
    end

    return solglouton, dot(solglouton,couts)
end
