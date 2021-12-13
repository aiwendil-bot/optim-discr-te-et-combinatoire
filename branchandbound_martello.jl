using LinearAlgebra

function backtracking(objets)
    variables_fixees = findall(x -> x == 1, objets)
    return variables_fixees[end]

end

function branchandbound(couts::Vector{Float64},poids::Vector{Float64}, capacite::Int64)
    compteurnoeuds::Int64 = 0
    bornemax::Int64 = bornemax = calcul_borne_martello(couts, poids,capacite)
    println("borne max : ", bornemax)
    bornemin = 0
    solrelax = zeros(length(poids))
    solglouton = zeros(length(poids))
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
        println(soltemp)
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

        println("borne : ",borne_LPK)

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
    return solglouton, dot(solglouton,couts), compteurnoeuds
end

function calcul_borne_martello(couts::Vector{Float64},poids::Vector{Float64}, capacite::Int64)::Int64
    solrelax = zeros(length(poids))
    solglouton = zeros(length(poids))

    s::Int64 = 1
    sommepoids::Int64 = 0
    while sommepoids + poids[s] <= capacite && s <= length(poids)
        sommepoids += poids[s]
        s += 1
        if s == length(poids) + 1
            return sum(couts)
        end
    end
    println(s)
    c_barre = capacite - sum([poids[i] for i in 1:(s-1)])
    U_0 = sum([couts[i] for i in 1:(s-1)]) + floor(c_barre * couts[s + 1]/poids[s + 1])
    U_1 = sum([couts[i] for i in 1:(s-1)]) + floor(couts[s] - (poids[s] - c_barre) * couts[s - 1] / poids[s - 1] )
    return max(U_0, U_1)
end
#instance 3
#println(branchandbound([4.0,9.0,10.0,9.0,3.0,14.0,14.0,2.0],[1.0,3.0,4.0,4.0,2.0,13.0,17.0,3.0],22))
#instance 1
#println(branchandbound([4.0,9.0,10.0,9.0,3.0,2.0],[1.0,3.0,4.0,4.0,2.0,3.0],7))
#instance 2
#println(branchandbound([70.0,20.0,39.0,37.0,7.0,5.0,10.0],[31.0,10.0,20.0,19.0,4.0,3.0,6.0],50))
#instance 4
println(branchandbound([112.0,90.0,15.0,12.0,12.0,9.0,26.0],[16.0,15.0,3.0,3.0,4.0,3.0,13.0],35))
