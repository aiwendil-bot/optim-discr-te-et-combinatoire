#=
HungFisk_BranchandBound:
- Julia version: 1.6.0
- Author: adrien callico & nicolas compère
- Date: 2021-12-04
=#

#=
pseudo code, notations :

F = ensemble des objets non assignés ie. quand F est vide, on a une solution admissible
S = ensemble des objets qui ont été assignés

STEP 1 (initialisation)

z = - infini, S vide, F = tous les objets et k=1

STEP 2 : (évaluation des bornes)

résoudre le problème relaxé pour les objets non assignés
si sol <= z* backtracking
si faisable (tous objets assignés) mettre à jour z* puis backtracking
sinon branching

STEP 3 (branching)

Sélectionner un objet i et l attribuer à un des sacs (y compris dummy)
évaluer la valeur des objets dans les sacs (sans le dummy)
si tous les objets ont été attribués, mettre à jour z* puis backtracking
sinon, mettre à jour S, F et k et recommencer cette étape

STEP 4 (mettre à jour z*) ok

STEP 5 (backtracking)

trouver le plus petit k tel que z_k0 (sol pb relaxé) <= z*
denoter le lvl précédent k_0 comme k_-1,
et dénoter l'objet correspondant dans S comme i_{k-1}'
si k_-1 <= 0, STOP
sinon
k = k_-1
et rendre libre tous les objets suivant i_{k-1}, et assigner l'objet i_{k-1} à un sac différent
et ré-évaluer (step 2)
s'il a déjà été attribué aux n+1 sacs, z = - inf et backtracking

précisions :

l'objet qu'on choisit d'assigner est celui de plus petit indice
les objets sont groupés par classe = de mêmes poids
les sacs sont groupés par classe = de mêmes capacités

un objet i peut être assigné à un sac j si soit
l'objet i est le premier de sa classe à être considéré
ou
le précédent objet a été assigné à un sac d'indice <= j

un objet i peut être assigné à un sac j si soit
j est le plus petit indice dans sa classe
ou
le sac j-1 n'est pas vide

un objet i peut être assigné au sac j si w_i <= capacité résiduelle du sac j

on stocke dans un vecteur u_i les sacs où un objet i peut être attribués
on prend le sac de plus petit indice
=#

include("solver_M01KP.jl")
include("generateInstance.jl")
using LinearAlgebra

function branchandbound_HungFisk(couts::Vector{Float64},poids::Vector{Float64}, capacites::Vector{Int64})
    nb_sacs = length(capacites)
    nb_objets = length(couts)
    capacites_residuelles = deepcopy(capacites)
    res = [zeros(Float64, nb_objets) for i in 1:(nb_sacs + 1)]
    coeff_optimal::Float64 = calcul_coeff_optimal(nb_objets, couts, poids, capacites) #relation 19 page 3

    #step 1
    res_etoile = [zeros(Int, nb_objets) for i in 1:(nb_sacs + 1)] #dummy knapsack
    z_etoile::Float64 = -Inf
    S = Vector{Int64}(undef,0)
    F = [i for i in 1:length(couts)]
    k::Int64 = 1
    tableau_backtracking = Vector{Float64}(undef,0) #valeurs = bornes duales
    bounding::Bool, branching::Bool, backtracking::Bool  = true, false, false

    while true

        #step 2
        if bounding
            bounding = false
            sol_relax, z_relax = solve_modelM01KP_surrogate(nb_objets, nb_sacs, couts, poids, capacites, coeff_optimal, S)
            push!(tableau_backtracking, z_relax)
            if z_relax <= z_etoile
                backtracking = true
                continue
            end
            feasible::Bool = true

            for j in 1:nb_sacs
                for l in 1:nb_objets
                    if sol_relax[j][l] != 0.0 || sol_relax[j][l] != 1.0
                        feasible = false
                        break
                    end
                end
            end
            #step 4
            if feasible
                z_etoile = z_relax
                res_etoile = sol_relax
                backtracking = true
            else
                branching = true
            end
            continue
        end

        #step 3
        if branching
            branching = false
            u = [Vector{Int64}(undef,0) for i in 1:nb_objets] #stocke les sacs à dos possibles pour chaque objet
            while length(F) > 0
                i = F[1] #on prend l'objet de plus petit indice
                u = calcul_sacs_disponibles(i, nb_sacs, capacites, poids, res, u)
                res[u[i][1]][i] = 1.0
                deleteat!(u[i],1)
                deleteat!(F, findall(x->x==i, F))
                push!(S,i)
                k += 1
            end
            #step 4
            if sum([dot(res[i],couts) for i in 1:nb_sacs]) > z_etoile
                z_etoile = sum([dot(res[i],couts) for i in 1:nb_sacs])
                res_etoile = copy(res)
                backtracking =true
            else
                branching = true
            end
            continue
        end
        #step 5
        if backtracking
            backtracking = false
            k = backtrack(tableau_backtracking, z_etoile)
            if k <= 0
                return res_etoile[1:nb_sacs], z_etoile
            else
                for l in (k+1):length(nb_objets)
                    deleteat!(S, findall(x->x==l, S))
                    push!(F, l)
                end
                sort!(F)
                if length(u[k]) > 0
                    for idx in 1:nb_sacs
                        if res[idx][k] == 1.0
                            res[idx][k] = 0.0
                            break
                        end
                    end
                    res[u[k][1]][k] = 1.0
                    bounding = true
                else
                    tableau_backtracking[k] = -Inf
                    backtracking = true
                end
                continue
            end
        end
    end
end

function backtrack(tableau_backtracking::Vector{Float64}, z_etoile::Float64)::Int64
    k_0::Int64 = 0
    for k in 1:length(tableau_backtracking)
        if tableau_backtracking[k] <= z_etoile
            k_0 = k
            break
        end
    end
    k_prec::Int64 = k_0 - 1
    return k_prec
end

function calcul_sacs_disponibles(i::Int64,nb_sacs::Int64,capacites::Vector{Int64},poids::Vector{Float64},res,u::Vector{Vector{Int64}})::Vector{Vector{Int64}}
    classe_objet = findall(x -> x==poids[i],poids)
    indice_dans_la_classe = findfirst(x->x==i,classe_objet)

    for j in 1:nb_sacs
        classe_sac = findall(x -> x==capacites[j], capacites)
        f_j = capacites[j] - dot(res[j],poids)
        attribue_not_greater::Bool = false
        for k in 1:j
            if res[k][i] == 1.0
                attribue_not_greater = true
            end
        end
        if (i == classe_objet[1] || !attribue_not_greater) && #rule 4
            (j == classe_sac[1] || dot(res[j],poids) > 0 ) && #=rule 5 =# (poids[i] <= f_j) #rule 6
            push!(u[i], j)
        end
    end
    if length(u[i]) == 0
        push!(u[i], nb_sacs+1)
    end
    return u
end

function calcul_coeff_optimal(nb_objets::Int64, couts::Vector{Float64}, poids::Vector{Float64}, capacites::Vector{Int64})::Float64
    somme_capacites = sum(capacites)
    somme_poids::Int64 = 0
    i::Int64 = 1
    if sum(poids) <= sum(capacites)
        return couts[end]/poids[end]
    end

    while somme_poids + poids[i] <= somme_capacites
        somme_poids += poids[i]
        i+= 1
    end
    return couts[i]/poids[i]
end
