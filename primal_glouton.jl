#= CALLICO / COMPÈRE - Université de Nantes - ORO - Optimisation discrète et combinatoire
Calcul d'une borne primale à l'aide d'un algorithme glouton pour le problème M01KP =#

using LinearAlgebra


function greedy(couts::Vector{Float64},poids::Vector{Float64}, capacites::Vector{Int64})
    nb_objets = length(couts)
    nb_sacs = length(capacites)
    capacites_residuelles = copy(capacites)
    res = [zeros(Int, nb_objets) for i in 1:nb_sacs]

    i, j = 1, 1

    while(i <= nb_sacs && j <= nb_objets)
        if poids[j] <= capacites_residuelles[i]
            res[i][j] = 1
            capacites_residuelles[i] -= poids[j]
            j += 1
        else
            j += 1
            i += 1
        end
    end
    return (res, dot(couts,res))
end
