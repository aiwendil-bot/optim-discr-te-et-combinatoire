#= CALLICO / COMPÈRE - Université de Nantes - ORO - Optimisation discrète et combinatoire
Calcul d'une borne primal à l'aide d'un algorithme glouton pour le problème M01KP =#

using Statistics

function item_critique(n::Int64, capa::Int64, couts::Vector{Int64}, poids::Vector{Int64})
    J1::Vector{Int64}
    J0::Vector{Int64}
    JC::Vector{Int64} = [1:n]
    capa_res::Int64 = capa
    partition::Bool = false
    while partition = false
        R::Vector{Int64} = [couts[j] / poids[j] for j in JC]
        λ::Float64 = median(R)
        G::Vector{Int64}
        L::Vector{Int64}
        E::Vector{Int64}
        for j in JC
            if couts[j] / poids[j] > λ
                push!(G, j)
            elseif couts[j] / poids[j] < λ
                push!(L, j)
            else
                push!(E, j)
            end
            c_p::Int64 = sum(poids[j] for j in G)
            c_pp::Int64 = c_p + sum(poids[j] for j in E)
            if c_p <= capa_res <= c_pp
                partition = true
            elseif c_p > capa_res
                J0 = J0
            end
        end
    end
end