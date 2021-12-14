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
include("Vb_relax_surrogate.jl")
using LinearAlgebra
using StatsBase
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
    tableau_backtracking = Dict{Int64,Float64}() #valeurs = bornes duales
    bounding::Bool, branching::Bool, backtracking::Bool  = true, false, false
    u = Vector{Vector{Int64}}(undef, 0)
    assignations = zeros(Int64,length(couts))
    deja_assignes = [Int64[] for i in 1:nb_objets]
    #z_relax = borneduale_surrogate(nb_sacs, assignations,couts, poids, capa, coeff, F)
    u = [Vector{Int64}(undef,0) for i in 1:nb_objets] #stocke les sacs à dos possibles pour chaque objet
    from_backtracking::Bool = false
    while true
        #step 2
        if bounding
            bounding = false
            #= essai d'ajout des règles 1 2

            couts_surrogate = Vector{Float64}(undef,0)
            poids_surrogate = Vector{Float64}(undef,0)
            capacites_surrogate = deepcopy(capacites)
            F_surrogate = Vector{Int64}(undef,0)
            F_calculs = deepcopy(F)
            min_poids_objets_libres = minimum([poids[i] for i in F])
            max_fj::Int64 = 0
            for j in 1:nb_sacs
                f_j = capacites[j] - dot(res[j],poids)
                if f_j < min_poids_objets_libres
                    deleteat!(capacites_surrogate,j)
                end
                if f_j > max_fj
                    max_fj = f_j
                end
            end
            for i in 1:length(F_calculs)
                if poids[F_calculs[i]] <= max_fj && length(F_surrogate) >0
                    push!(couts_surrogate,couts[F_calculs[i]])
                    push!(poids_surrogate,poids[F_calculs[i]])
                    for j in (i+1):length(F_calculs)
                        F_calculs[j] -= 1
                    end
                    deleteat!(F_surrogate,findfirst(x->x==i,F_calculs))
                end
            end
            z_relax = borneduale_surrogate(couts_surrogate, poids_surrogate, capacites_surrogate, coeff_optimal, F_surrogate)
            fin du test d'ajout des règles =#




            z_relax = borneduale_surrogate(length(capacites), assignations,couts, poids, capacites, coeff_optimal, F)
            tableau_backtracking[k] = z_relax
            if z_relax <= z_etoile
                backtracking = true
                continue
            end
            #step 4
            branching = true
            continue
        end

        #step 3
        if branching
            branching = false

            while length(F) > 0

                i = F[1] #on prend l'objet de plus petit indice
                u = calcul_sacs_disponibles(deja_assignes,i, nb_sacs, capacites, poids, res, u)

                res[u[i][1]][i] = 1.0
                assignations[i] = u[i][1]
                push!(deja_assignes[i],u[i][1])
                deleteat!(u[i],1)
                deleteat!(F, findall(x->x==i, F))
                push!(S,i)
                k += 1

                z_relax = borneduale_surrogate(length(capacites),assignations,couts, poids, capacites, coeff_optimal, F)
                tableau_backtracking[k] = z_relax
                #sol_relax, z_relax = solve_modelM01KP_surrogate(nb_objets, nb_sacs, couts, poids, capacites, coeff_optimal, S)
                #z_relax = borneduale_surrogate(couts, poids, capacites, coeff_optimal, S)
                #tableau_backtracking[k] = z_relax
                 #       println(S)
                #println(tableau_backtracking)
            end
            #step 4
            if sum([dot(res[i],couts) for i in 1:nb_sacs]) > z_etoile
                z_etoile = sum([dot(res[i],couts) for i in 1:nb_sacs])
                res_etoile = copy(res)

            end

            backtracking =true
            continue
        end
        #step 5
        if backtracking
            backtracking = false
            k_prec = k
            k = backtrack(tableau_backtracking, z_etoile)
            if k <= 0
                return res_etoile[1:nb_sacs], z_etoile
            else
                for l in (k+1):nb_objets
                    deleteat!(S, findall(x->x==l, S))
                    push!(F, l)
                    assignations[l]=0

                end
                F = [key for (key, val) in countmap(F) if val >= 1]
                sort!(F)
                for idx in 1:nb_sacs, j in (k+1):nb_objets
                        if res[idx][j] == 1.0
                            res[idx][j] = 0.0
                        end
                end
                if length(u[k]) > 0

                    for i in 1:length(capacites)
                        res[i][k] = 0
                    end
                    res[u[k][1]][k] = 1.0
                    assignations[k] = u[k][1]
                    push!(deja_assignes[k],u[k][1])
                    deleteat!(u[k],1)
                    push!(S,k)
                    deleteat!(F, findall(x->x==k, F))
                    bounding = true
                    from_backtracking = true

                else
                    if length(findall(x->x==(nb_sacs+1),deja_assignes[k])) == 0
                        assignations[k] = nb_sacs + 1
                        push!(deja_assignes[k],nb_sacs + 1)
                        push!(S,k)
                        deleteat!(F, findall(x->x==k, F))
                        for i in 1:length(capacites)
                            res[i][k] = 0
                        end
                        bounding = true
                        from_backtracking = true

                    else
                        tableau_backtracking[k] = - 2^63

                        backtracking = true

                    end
                end

                continue
            end
        end
    end
end

function backtrack(tableau_backtracking::Dict{Int64,Float64}, z_etoile::Float64)::Int64
    k_0::Int64 = 2^63 - 1
    min_k = k_0
    for (k,v) in tableau_backtracking
        if v <= z_etoile
            if k < min_k
                min_k = k
            end
        end
    end
    return min_k - 1
end

function calcul_sacs_disponibles(deja_assignes::Vector{Vector{Int64}}, i::Int64,nb_sacs::Int64,capacites::Vector{Int64},poids::Vector{Float64},res,u::Vector{Vector{Int64}})::Vector{Vector{Int64}}
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

        if (poids[i] <= f_j) && length(findall(x->x==j,deja_assignes[i])) == 0  #(i == classe_objet[1] || !attribue_not_greater) && #rule 4
            #=(j == classe_sac[1] || dot(res[j],poids) > 0 ) && rule 5 =#  #rule 6
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

#non utilisée
function isFeasible(n::Int64, m::Int64, sol::Matrix{Float64}, poids::Vector{Float64}, capacites::Vector{Int64})::Bool
    for i in 1:m
        if sum([sol[i,j] * poids[j] for j in 1:n]) > capacites[i]
            return false
        end
    end
    return true
end

instance1 = instance(10, 2, [44.0, 51.0, 37.0, 97.0, 57.0, 86.0, 96.0, 40.0, 24.0, 1.0], [2.0, 3.0, 10.0, 59.0, 39.0, 60.0, 84.0, 52.0, 87.0, 47.0], [27, 19])
instance2 = instance(10, 7, [95.0, 43.0, 49.0, 100.0, 33.0, 17.0, 29.0, 33.0, 42.0, 12.0], [23.0, 16.0, 22.0, 77.0, 42.0, 23.0, 42.0, 60.0, 96.0, 61.0], [28, 22, 18, 14, 13, 13, 13])
instance3 = instance(50, 5, [100.0, 34.0, 75.0, 53.0, 62.0, 93.0, 35.0, 75.0, 77.0, 94.0, 85.0, 80.0, 67.0, 76.0, 41.0, 34.0, 90.0, 38.0, 56.0, 91.0, 87.0, 48.0, 75.0, 80.0, 87.0, 84.0, 59.0, 84.0, 36.0, 18.0, 61.0, 39.0, 77.0, 77.0, 33.0, 46.0, 68.0, 26.0, 22.0, 13.0, 50.0, 38.0, 34.0, 35.0, 39.0, 13.0, 12.0, 16.0, 3.0, 2.0], [2.0, 5.0, 12.0, 9.0, 16.0, 26.0, 10.0, 22.0, 23.0, 30.0, 30.0, 34.0, 29.0, 35.0, 20.0, 19.0, 52.0, 23.0, 35.0, 67.0, 67.0, 38.0, 66.0, 72.0, 79.0, 86.0, 61.0, 88.0, 38.0, 20.0, 69.0, 45.0, 93.0, 98.0, 43.0, 64.0, 96.0, 38.0, 37.0, 22.0, 95.0, 74.0, 79.0, 83.0, 96.0, 32.0, 53.0, 85.0, 25.0, 73.0], [22, 20, 13, 11, 11])
instance4 = instance(50, 30, [57.0, 86.0, 50.0, 48.0, 70.0, 55.0, 54.0, 75.0, 69.0, 68.0, 78.0, 83.0, 97.0, 63.0, 66.0, 76.0, 58.0, 51.0, 99.0, 100.0, 93.0, 60.0, 87.0, 51.0, 67.0, 87.0, 63.0, 87.0, 43.0, 72.0, 26.0, 66.0, 70.0, 65.0, 54.0, 58.0, 35.0, 12.0, 24.0, 44.0, 14.0, 36.0, 31.0, 14.0, 18.0, 1.0, 7.0, 4.0, 2.0, 1.0], [5.0, 15.0, 14.0, 14.0, 21.0, 22.0, 22.0, 31.0, 29.0, 33.0, 42.0, 50.0, 64.0, 42.0, 48.0, 61.0, 48.0, 45.0, 90.0, 91.0, 85.0, 57.0, 84.0, 50.0, 67.0, 87.0, 64.0, 93.0, 48.0, 83.0, 30.0, 82.0, 87.0, 96.0, 91.0, 98.0, 64.0, 22.0, 50.0, 98.0, 35.0, 100.0, 93.0, 57.0, 85.0, 10.0, 76.0, 63.0, 41.0, 24.0], [30, 28, 28, 28, 27, 27, 27, 25, 25, 24, 24, 23, 22, 22, 21, 20, 20, 17, 17, 17, 14, 14, 14, 14, 13, 12, 11, 10, 10, 10])
instance5 = instance(100, 25, [73.0, 70.0, 44.0, 85.0, 79.0, 33.0, 86.0, 97.0, 84.0, 90.0, 92.0, 67.0, 78.0, 64.0, 31.0, 76.0, 38.0, 86.0, 60.0, 85.0, 78.0, 95.0, 55.0, 19.0, 83.0, 87.0, 33.0, 100.0, 98.0, 32.0, 84.0, 98.0, 61.0, 50.0, 23.0, 56.0, 62.0, 90.0, 96.0, 90.0, 50.0, 83.0, 73.0, 59.0, 96.0, 92.0, 90.0, 53.0, 75.0, 71.0, 28.0, 92.0, 49.0, 25.0, 13.0, 15.0, 55.0, 96.0, 86.0, 55.0, 14.0, 82.0, 90.0, 87.0, 32.0, 29.0, 45.0, 16.0, 79.0, 76.0, 61.0, 47.0, 72.0, 38.0, 59.0, 27.0, 45.0, 38.0, 20.0, 3.0, 40.0, 25.0, 18.0, 35.0, 26.0, 35.0, 15.0, 14.0, 25.0, 26.0, 16.0, 14.0, 18.0, 14.0, 14.0, 8.0, 3.0, 3.0, 3.0, 2.0], [3.0, 3.0, 2.0, 5.0, 5.0, 3.0, 11.0, 15.0, 13.0, 15.0, 19.0, 14.0, 17.0, 14.0, 7.0, 20.0, 11.0, 25.0, 18.0, 34.0, 32.0, 39.0, 23.0, 8.0, 36.0, 41.0, 17.0, 53.0, 52.0, 17.0, 46.0, 55.0, 35.0, 29.0, 14.0, 35.0, 39.0, 57.0, 66.0, 63.0, 36.0, 64.0, 57.0, 49.0, 80.0, 77.0, 76.0, 45.0, 66.0, 63.0, 25.0, 83.0, 45.0, 23.0, 12.0, 14.0, 53.0, 93.0, 87.0, 56.0, 15.0, 88.0, 97.0, 97.0, 36.0, 33.0, 53.0, 19.0, 96.0, 93.0, 75.0, 62.0, 95.0, 57.0, 100.0, 50.0, 85.0, 74.0, 40.0, 6.0, 92.0, 58.0, 42.0, 83.0, 66.0, 93.0, 41.0, 45.0, 88.0, 94.0, 72.0, 64.0, 89.0, 74.0, 98.0, 79.0, 49.0, 64.0, 83.0, 100.0], [30, 30, 29, 27, 27, 27, 26, 25, 23, 22, 21, 21, 19, 19, 19, 18, 17, 16, 13, 12, 12, 12, 12, 11, 10])
instance6 = instance(100, 75, [47.0, 65.0, 46.0, 75.0, 45.0, 54.0, 34.0, 29.0, 83.0, 83.0, 76.0, 91.0, 93.0, 77.0, 36.0, 75.0, 89.0, 93.0, 27.0, 81.0, 51.0, 67.0, 79.0, 98.0, 78.0, 88.0, 95.0, 70.0, 56.0, 89.0, 52.0, 63.0, 92.0, 57.0, 59.0, 52.0, 92.0, 56.0, 94.0, 97.0, 64.0, 73.0, 92.0, 43.0, 89.0, 98.0, 54.0, 42.0, 97.0, 33.0, 83.0, 62.0, 56.0, 91.0, 65.0, 59.0, 34.0, 81.0, 69.0, 79.0, 24.0, 56.0, 52.0, 44.0, 75.0, 51.0, 68.0, 15.0, 52.0, 64.0, 26.0, 36.0, 10.0, 22.0, 53.0, 38.0, 51.0, 40.0, 21.0, 31.0, 35.0, 21.0, 34.0, 33.0, 8.0, 30.0, 25.0, 12.0, 22.0, 17.0, 9.0, 15.0, 20.0, 19.0, 3.0, 18.0, 15.0, 2.0, 5.0, 2.0], [1.0, 3.0, 7.0, 15.0, 9.0, 12.0, 8.0, 7.0, 24.0, 26.0, 24.0, 30.0, 32.0, 27.0, 14.0, 30.0, 36.0, 38.0, 12.0, 38.0, 27.0, 36.0, 43.0, 55.0, 45.0, 53.0, 58.0, 43.0, 35.0, 56.0, 33.0, 41.0, 62.0, 42.0, 46.0, 42.0, 77.0, 47.0, 79.0, 82.0, 56.0, 66.0, 84.0, 40.0, 84.0, 96.0, 53.0, 42.0, 98.0, 34.0, 87.0, 65.0, 60.0, 100.0, 76.0, 69.0, 40.0, 96.0, 82.0, 95.0, 30.0, 71.0, 67.0, 58.0, 100.0, 71.0, 95.0, 21.0, 76.0, 100.0, 42.0, 61.0, 17.0, 39.0, 96.0, 73.0, 99.0, 85.0, 47.0, 70.0, 80.0, 51.0, 86.0, 85.0, 21.0, 79.0, 78.0, 43.0, 85.0, 66.0, 37.0, 63.0, 85.0, 86.0, 14.0, 97.0, 81.0, 19.0, 51.0, 75.0], [30, 30, 30, 29, 29, 29, 29, 29, 29, 28, 28, 27, 27, 27, 27, 26, 26, 26, 26, 26, 26, 26, 25, 25, 25, 25, 24, 24, 24, 24, 22, 22, 22, 22, 21, 21, 21, 20, 19, 19, 19, 19, 19, 18, 18, 18, 18, 18, 18, 18, 17, 17, 16, 16, 16, 15, 14, 14, 14, 14, 14, 13, 13, 13, 13, 13, 13, 12, 11, 10, 10, 10, 10, 10, 10])
instance7 = instance(200, 50, [94.0, 59.0, 88.0, 22.0, 68.0, 47.0, 44.0, 91.0, 99.0, 58.0, 68.0, 38.0, 69.0, 70.0, 34.0, 93.0, 93.0, 82.0, 55.0, 48.0, 75.0, 68.0, 99.0, 46.0, 77.0, 69.0, 90.0, 90.0, 72.0, 11.0, 83.0, 67.0, 77.0, 71.0, 42.0, 84.0, 71.0, 58.0, 44.0, 46.0, 83.0, 33.0, 90.0, 83.0, 70.0, 64.0, 71.0, 84.0, 87.0, 32.0, 69.0, 73.0, 63.0, 64.0, 60.0, 99.0, 94.0, 50.0, 70.0, 62.0, 92.0, 54.0, 36.0, 61.0, 22.0, 80.0, 12.0, 86.0, 87.0, 46.0, 86.0, 95.0, 95.0, 20.0, 52.0, 95.0, 44.0, 82.0, 47.0, 41.0, 77.0, 83.0, 74.0, 89.0, 89.0, 46.0, 84.0, 85.0, 88.0, 97.0, 52.0, 48.0, 92.0, 86.0, 79.0, 51.0, 34.0, 54.0, 74.0, 53.0, 52.0, 85.0, 61.0, 70.0, 69.0, 47.0, 72.0, 64.0, 45.0, 53.0, 53.0, 56.0, 60.0, 60.0, 15.0, 45.0, 58.0, 50.0, 41.0, 22.0, 8.0, 55.0, 26.0, 21.0, 36.0, 53.0, 33.0, 44.0, 44.0, 53.0, 44.0, 35.0, 27.0, 35.0, 41.0, 45.0, 35.0, 21.0, 32.0, 40.0, 33.0, 30.0, 40.0, 27.0, 27.0, 28.0, 28.0, 24.0, 36.0, 9.0, 29.0, 25.0, 31.0, 19.0, 18.0, 26.0, 26.0, 28.0, 28.0, 26.0, 26.0, 24.0, 15.0, 4.0, 12.0, 17.0, 11.0, 24.0, 6.0, 5.0, 18.0, 18.0, 22.0, 14.0, 12.0, 18.0, 10.0, 7.0, 8.0, 14.0, 9.0, 3.0, 2.0, 7.0, 8.0, 11.0, 3.0, 7.0, 8.0, 4.0, 8.0, 7.0, 6.0, 4.0, 5.0, 3.0, 2.0, 3.0, 2.0, 1.0], [1.0, 2.0, 3.0, 1.0, 4.0, 4.0, 4.0, 9.0, 13.0, 8.0, 10.0, 7.0, 15.0, 16.0, 8.0, 22.0, 23.0, 21.0, 16.0, 14.0, 22.0, 20.0, 30.0, 14.0, 24.0, 22.0, 29.0, 29.0, 24.0, 4.0, 32.0, 26.0, 30.0, 28.0, 17.0, 35.0, 30.0, 25.0, 19.0, 21.0, 38.0, 16.0, 44.0, 41.0, 35.0, 32.0, 36.0, 44.0, 46.0, 17.0, 37.0, 40.0, 36.0, 38.0, 36.0, 60.0, 58.0, 31.0, 44.0, 39.0, 58.0, 36.0, 24.0, 43.0, 16.0, 60.0, 9.0, 65.0, 66.0, 35.0, 66.0, 75.0, 75.0, 16.0, 42.0, 78.0, 38.0, 72.0, 42.0, 37.0, 71.0, 77.0, 69.0, 84.0, 84.0, 44.0, 81.0, 82.0, 85.0, 98.0, 53.0, 49.0, 94.0, 91.0, 86.0, 57.0, 38.0, 61.0, 86.0, 62.0, 61.0, 100.0, 74.0, 86.0, 85.0, 58.0, 90.0, 81.0, 58.0, 70.0, 73.0, 81.0, 87.0, 90.0, 23.0, 69.0, 89.0, 77.0, 64.0, 35.0, 13.0, 90.0, 43.0, 35.0, 60.0, 89.0, 56.0, 75.0, 75.0, 92.0, 80.0, 64.0, 51.0, 68.0, 80.0, 92.0, 75.0, 45.0, 69.0, 87.0, 72.0, 69.0, 94.0, 64.0, 67.0, 70.0, 70.0, 61.0, 99.0, 25.0, 81.0, 70.0, 87.0, 59.0, 58.0, 84.0, 86.0, 95.0, 99.0, 94.0, 96.0, 89.0, 56.0, 15.0, 46.0, 68.0, 44.0, 99.0, 25.0, 21.0, 79.0, 81.0, 100.0, 65.0, 57.0, 88.0, 51.0, 37.0, 43.0, 77.0, 52.0, 21.0, 15.0, 53.0, 65.0, 96.0, 27.0, 68.0, 78.0, 41.0, 93.0, 95.0, 86.0, 60.0, 79.0, 59.0, 53.0, 94.0, 64.0, 82.0], [30, 30, 30, 29, 29, 29, 29, 28, 28, 28, 27, 26, 26, 26, 26, 26, 25, 25, 25, 23, 22, 22, 21, 21, 20, 20, 20, 19, 18, 18, 18, 18, 18, 17, 16, 16, 14, 14, 14, 14, 13, 13, 12, 12, 12, 12, 12, 10, 10, 10])
instance8 = instance(200, 125, [53.0, 75.0, 61.0, 83.0, 53.0, 90.0, 98.0, 69.0, 89.0, 82.0, 54.0, 70.0, 57.0, 16.0, 79.0, 75.0, 22.0, 80.0, 61.0, 64.0, 92.0, 84.0, 97.0, 61.0, 38.0, 76.0, 60.0, 64.0, 21.0, 97.0, 88.0, 94.0, 100.0, 56.0, 91.0, 84.0, 93.0, 90.0, 32.0, 75.0, 90.0, 88.0, 57.0, 36.0, 75.0, 91.0, 85.0, 68.0, 64.0, 99.0, 67.0, 65.0, 93.0, 27.0, 52.0, 95.0, 73.0, 26.0, 53.0, 95.0, 81.0, 71.0, 81.0, 77.0, 97.0, 40.0, 81.0, 99.0, 45.0, 88.0, 83.0, 31.0, 81.0, 81.0, 93.0, 60.0, 45.0, 95.0, 77.0, 94.0, 97.0, 99.0, 75.0, 55.0, 91.0, 62.0, 97.0, 85.0, 70.0, 33.0, 66.0, 60.0, 44.0, 97.0, 15.0, 33.0, 84.0, 85.0, 78.0, 94.0, 75.0, 83.0, 29.0, 26.0, 70.0, 91.0, 52.0, 83.0, 12.0, 81.0, 75.0, 75.0, 70.0, 55.0, 82.0, 52.0, 81.0, 69.0, 76.0, 50.0, 75.0, 12.0, 24.0, 17.0, 25.0, 70.0, 55.0, 50.0, 19.0, 61.0, 9.0, 45.0, 68.0, 41.0, 50.0, 59.0, 22.0, 22.0, 28.0, 9.0, 41.0, 14.0, 58.0, 38.0, 33.0, 55.0, 35.0, 45.0, 53.0, 39.0, 18.0, 43.0, 37.0, 33.0, 43.0, 23.0, 33.0, 34.0, 27.0, 17.0, 35.0, 21.0, 15.0, 31.0, 35.0, 26.0, 30.0, 26.0, 14.0, 31.0, 22.0, 12.0, 18.0, 25.0, 23.0, 11.0, 20.0, 10.0, 4.0, 18.0, 18.0, 17.0, 6.0, 14.0, 15.0, 9.0, 9.0, 6.0, 9.0, 3.0, 8.0, 4.0, 8.0, 8.0, 7.0, 6.0, 3.0, 6.0, 5.0, 2.0], [1.0, 2.0, 2.0, 3.0, 2.0, 4.0, 5.0, 4.0, 6.0, 7.0, 5.0, 8.0, 7.0, 2.0, 10.0, 10.0, 3.0, 11.0, 9.0, 13.0, 19.0, 19.0, 22.0, 14.0, 9.0, 18.0, 15.0, 17.0, 6.0, 28.0, 26.0, 28.0, 30.0, 17.0, 28.0, 27.0, 31.0, 30.0, 11.0, 27.0, 35.0, 35.0, 23.0, 15.0, 32.0, 40.0, 38.0, 32.0, 31.0, 48.0, 34.0, 33.0, 48.0, 15.0, 29.0, 53.0, 42.0, 15.0, 31.0, 56.0, 48.0, 43.0, 50.0, 49.0, 63.0, 26.0, 53.0, 65.0, 30.0, 59.0, 56.0, 21.0, 56.0, 58.0, 67.0, 46.0, 35.0, 74.0, 60.0, 74.0, 77.0, 80.0, 61.0, 45.0, 76.0, 52.0, 82.0, 74.0, 61.0, 29.0, 59.0, 54.0, 40.0, 90.0, 14.0, 31.0, 79.0, 80.0, 74.0, 90.0, 72.0, 80.0, 29.0, 26.0, 73.0, 96.0, 55.0, 89.0, 13.0, 90.0, 84.0, 85.0, 81.0, 64.0, 97.0, 62.0, 98.0, 87.0, 97.0, 65.0, 100.0, 16.0, 32.0, 23.0, 34.0, 97.0, 77.0, 70.0, 27.0, 87.0, 13.0, 65.0, 99.0, 60.0, 76.0, 90.0, 36.0, 36.0, 46.0, 15.0, 70.0, 24.0, 100.0, 67.0, 60.0, 100.0, 65.0, 84.0, 100.0, 74.0, 36.0, 88.0, 79.0, 74.0, 97.0, 53.0, 79.0, 82.0, 66.0, 42.0, 89.0, 56.0, 41.0, 85.0, 98.0, 75.0, 87.0, 76.0, 43.0, 100.0, 74.0, 41.0, 66.0, 96.0, 95.0, 51.0, 94.0, 48.0, 20.0, 93.0, 97.0, 94.0, 35.0, 82.0, 91.0, 76.0, 82.0, 56.0, 90.0, 31.0, 85.0, 47.0, 95.0, 96.0, 85.0, 87.0, 44.0, 95.0, 87.0, 47.0], [30, 30, 30, 30, 30, 30, 30, 30, 29, 29, 29, 29, 29, 29, 29, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 26, 26, 26, 26, 26, 26, 25, 25, 25, 25, 24, 24, 24, 24, 24, 24, 24, 23, 23, 23, 23, 23, 23, 22, 22, 22, 22, 22, 21, 21, 21, 21, 21, 21, 21, 21, 20, 20, 20, 20, 20, 20, 20, 19, 19, 19, 19, 19, 19, 18, 18, 18, 18, 18, 18, 17, 17, 16, 16, 16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 14, 14, 14, 14, 14, 14, 13, 13, 13, 13, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10])
instance9 = instance(100, 50, [87.0, 66.0, 57.0, 68.0, 56.0, 9.0, 95.0, 96.0, 70.0, 41.0, 95.0, 82.0, 40.0, 77.0, 54.0, 95.0, 99.0, 95.0, 98.0, 90.0, 55.0, 44.0, 58.0, 98.0, 91.0, 94.0, 74.0, 80.0, 54.0, 95.0, 85.0, 85.0, 84.0, 71.0, 69.0, 74.0, 24.0, 98.0, 72.0, 93.0, 78.0, 11.0, 42.0, 67.0, 75.0, 24.0, 51.0, 51.0, 46.0, 85.0, 31.0, 18.0, 88.0, 98.0, 76.0, 59.0, 52.0, 97.0, 35.0, 82.0, 69.0, 74.0, 89.0, 71.0, 38.0, 57.0, 50.0, 74.0, 70.0, 79.0, 35.0, 74.0, 60.0, 53.0, 23.0, 30.0, 54.0, 57.0, 54.0, 32.0, 50.0, 44.0, 40.0, 30.0, 34.0, 44.0, 38.0, 28.0, 28.0, 31.0, 25.0, 17.0, 12.0, 19.0, 8.0, 11.0, 4.0, 3.0, 1.0, 1.0], [1.0, 1.0, 1.0, 5.0, 6.0, 1.0, 11.0, 12.0, 10.0, 7.0, 18.0, 16.0, 8.0, 17.0, 13.0, 23.0, 24.0, 25.0, 26.0, 24.0, 16.0, 13.0, 18.0, 32.0, 31.0, 34.0, 27.0, 33.0, 24.0, 48.0, 45.0, 45.0, 45.0, 40.0, 39.0, 46.0, 15.0, 62.0, 47.0, 62.0, 53.0, 8.0, 31.0, 51.0, 58.0, 19.0, 41.0, 41.0, 40.0, 74.0, 27.0, 16.0, 80.0, 90.0, 72.0, 57.0, 52.0, 97.0, 35.0, 83.0, 72.0, 78.0, 97.0, 78.0, 43.0, 66.0, 59.0, 88.0, 84.0, 95.0, 45.0, 98.0, 82.0, 80.0, 35.0, 47.0, 92.0, 98.0, 97.0, 61.0, 98.0, 93.0, 86.0, 66.0, 76.0, 99.0, 89.0, 66.0, 67.0, 90.0, 86.0, 61.0, 45.0, 96.0, 59.0, 99.0, 60.0, 77.0, 75.0, 89.0], [30, 30, 30, 30, 29, 29, 28, 27, 27, 26, 26, 25, 25, 24, 24, 24, 23, 22, 22, 22, 22, 22, 20, 20, 20, 20, 19, 19, 18, 17, 17, 17, 17, 16, 15, 15, 15, 14, 14, 14, 14, 13, 13, 13, 12, 12, 11, 11, 11, 10])
instance10 = instance(150, 20, [80.0, 97.0, 47.0, 60.0, 77.0, 51.0, 43.0, 25.0, 61.0, 60.0, 39.0, 84.0, 41.0, 18.0, 56.0, 46.0, 70.0, 77.0, 88.0, 73.0, 42.0, 48.0, 48.0, 12.0, 50.0, 39.0, 11.0, 29.0, 42.0, 73.0, 46.0, 48.0, 59.0, 34.0, 69.0, 75.0, 68.0, 73.0, 72.0, 40.0, 20.0, 94.0, 93.0, 51.0, 46.0, 81.0, 49.0, 90.0, 92.0, 55.0, 63.0, 83.0, 80.0, 40.0, 91.0, 95.0, 86.0, 98.0, 56.0, 82.0, 67.0, 50.0, 59.0, 86.0, 91.0, 56.0, 25.0, 100.0, 67.0, 89.0, 97.0, 89.0, 54.0, 80.0, 48.0, 69.0, 80.0, 62.0, 90.0, 78.0, 47.0, 37.0, 51.0, 76.0, 51.0, 22.0, 79.0, 61.0, 67.0, 74.0, 77.0, 16.0, 19.0, 44.0, 46.0, 47.0, 39.0, 11.0, 48.0, 59.0, 54.0, 46.0, 14.0, 45.0, 43.0, 57.0, 39.0, 40.0, 29.0, 51.0, 28.0, 14.0, 31.0, 36.0, 34.0, 36.0, 25.0, 29.0, 43.0, 41.0, 30.0, 23.0, 15.0, 27.0, 28.0, 15.0, 30.0, 22.0, 7.0, 23.0, 20.0, 21.0, 15.0, 22.0, 20.0, 18.0, 5.0, 12.0, 10.0, 12.0, 14.0, 1.0, 4.0, 9.0, 6.0, 5.0, 1.0, 2.0, 1.0, 1.0], [2.0, 4.0, 2.0, 3.0, 4.0, 4.0, 5.0, 3.0, 8.0, 9.0, 6.0, 18.0, 9.0, 4.0, 13.0, 11.0, 18.0, 20.0, 23.0, 20.0, 12.0, 16.0, 16.0, 4.0, 17.0, 14.0, 4.0, 11.0, 16.0, 28.0, 19.0, 21.0, 29.0, 17.0, 38.0, 44.0, 40.0, 43.0, 43.0, 24.0, 12.0, 60.0, 61.0, 36.0, 33.0, 59.0, 36.0, 67.0, 69.0, 42.0, 50.0, 66.0, 64.0, 32.0, 73.0, 78.0, 73.0, 85.0, 49.0, 73.0, 60.0, 45.0, 54.0, 80.0, 86.0, 53.0, 24.0, 97.0, 65.0, 88.0, 99.0, 93.0, 57.0, 85.0, 52.0, 75.0, 87.0, 68.0, 99.0, 86.0, 52.0, 42.0, 59.0, 89.0, 60.0, 26.0, 98.0, 78.0, 86.0, 95.0, 99.0, 21.0, 25.0, 58.0, 61.0, 63.0, 53.0, 15.0, 67.0, 83.0, 77.0, 66.0, 21.0, 70.0, 67.0, 92.0, 64.0, 67.0, 49.0, 89.0, 49.0, 25.0, 56.0, 66.0, 66.0, 70.0, 50.0, 60.0, 92.0, 89.0, 66.0, 52.0, 37.0, 75.0, 85.0, 46.0, 95.0, 72.0, 26.0, 91.0, 80.0, 88.0, 65.0, 98.0, 90.0, 87.0, 28.0, 69.0, 58.0, 72.0, 92.0, 7.0, 32.0, 81.0, 67.0, 74.0, 31.0, 90.0, 50.0, 100.0], [29, 27, 25, 25, 24, 21, 21, 19, 19, 18, 18, 18, 16, 15, 15, 15, 15, 14, 10, 10])
instance11 = instance(150, 20, [84.0, 98.0, 60.0, 89.0, 77.0, 25.0, 58.0, 57.0, 33.0, 89.0, 47.0, 59.0, 71.0, 85.0, 99.0, 58.0, 66.0, 16.0, 85.0, 94.0, 77.0, 79.0, 95.0, 51.0, 19.0, 94.0, 28.0, 97.0, 93.0, 20.0, 69.0, 99.0, 42.0, 98.0, 98.0, 97.0, 39.0, 63.0, 56.0, 81.0, 76.0, 97.0, 88.0, 87.0, 93.0, 67.0, 28.0, 18.0, 44.0, 84.0, 100.0, 78.0, 75.0, 99.0, 92.0, 80.0, 85.0, 90.0, 96.0, 87.0, 63.0, 21.0, 100.0, 48.0, 58.0, 54.0, 8.0, 64.0, 55.0, 25.0, 27.0, 100.0, 44.0, 63.0, 98.0, 78.0, 68.0, 69.0, 68.0, 99.0, 91.0, 96.0, 66.0, 46.0, 89.0, 49.0, 61.0, 39.0, 67.0, 49.0, 26.0, 18.0, 35.0, 64.0, 34.0, 79.0, 65.0, 76.0, 50.0, 15.0, 21.0, 20.0, 63.0, 53.0, 58.0, 46.0, 58.0, 35.0, 43.0, 39.0, 46.0, 34.0, 43.0, 47.0, 32.0, 5.0, 38.0, 25.0, 7.0, 22.0, 19.0, 27.0, 25.0, 19.0, 15.0, 19.0, 12.0, 18.0, 13.0, 1.0, 14.0, 8.0, 14.0, 14.0, 15.0, 16.0, 16.0, 5.0, 9.0, 9.0, 2.0, 13.0, 3.0, 8.0, 9.0, 4.0, 2.0, 5.0, 3.0, 3.0], [1.0, 3.0, 2.0, 3.0, 3.0, 1.0, 3.0, 3.0, 2.0, 6.0, 5.0, 7.0, 9.0, 12.0, 15.0, 9.0, 11.0, 3.0, 16.0, 18.0, 15.0, 17.0, 22.0, 12.0, 5.0, 26.0, 8.0, 30.0, 29.0, 7.0, 25.0, 36.0, 17.0, 40.0, 40.0, 40.0, 17.0, 30.0, 27.0, 43.0, 41.0, 53.0, 51.0, 51.0, 55.0, 40.0, 17.0, 11.0, 27.0, 52.0, 64.0, 50.0, 49.0, 67.0, 63.0, 61.0, 66.0, 70.0, 75.0, 68.0, 50.0, 17.0, 82.0, 41.0, 50.0, 47.0, 7.0, 57.0, 50.0, 23.0, 25.0, 93.0, 41.0, 59.0, 93.0, 75.0, 66.0, 67.0, 67.0, 98.0, 91.0, 98.0, 68.0, 48.0, 94.0, 53.0, 67.0, 44.0, 76.0, 56.0, 30.0, 21.0, 41.0, 75.0, 40.0, 96.0, 79.0, 95.0, 64.0, 20.0, 28.0, 27.0, 86.0, 80.0, 89.0, 72.0, 94.0, 58.0, 75.0, 69.0, 82.0, 62.0, 83.0, 95.0, 69.0, 11.0, 84.0, 60.0, 19.0, 67.0, 58.0, 85.0, 79.0, 63.0, 50.0, 67.0, 43.0, 66.0, 51.0, 4.0, 58.0, 37.0, 66.0, 72.0, 78.0, 84.0, 88.0, 30.0, 56.0, 63.0, 14.0, 99.0, 28.0, 76.0, 100.0, 54.0, 27.0, 71.0, 44.0, 62.0], [28, 26, 24, 24, 24, 23, 22, 20, 19, 19, 19, 17, 17, 16, 16, 16, 15, 12, 12, 10])
instance12 = instance(300, 20, [96.0, 80.0, 67.0, 99.0, 31.0, 100.0, 53.0, 77.0, 75.0, 86.0, 54.0, 96.0, 93.0, 33.0, 96.0, 89.0, 48.0, 80.0, 85.0, 47.0, 62.0, 74.0, 59.0, 100.0, 78.0, 74.0, 22.0, 52.0, 85.0, 15.0, 76.0, 98.0, 76.0, 72.0, 96.0, 67.0, 78.0, 65.0, 29.0, 92.0, 44.0, 59.0, 89.0, 33.0, 18.0, 45.0, 68.0, 93.0, 20.0, 88.0, 64.0, 33.0, 91.0, 98.0, 80.0, 78.0, 29.0, 84.0, 44.0, 81.0, 64.0, 48.0, 35.0, 93.0, 58.0, 17.0, 36.0, 45.0, 69.0, 77.0, 79.0, 81.0, 85.0, 76.0, 92.0, 24.0, 59.0, 90.0, 88.0, 95.0, 21.0, 59.0, 90.0, 89.0, 100.0, 75.0, 82.0, 80.0, 88.0, 86.0, 55.0, 98.0, 33.0, 98.0, 81.0, 27.0, 55.0, 85.0, 68.0, 52.0, 63.0, 66.0, 96.0, 68.0, 88.0, 38.0, 6.0, 31.0, 87.0, 87.0, 97.0, 55.0, 83.0, 97.0, 41.0, 15.0, 34.0, 61.0, 61.0, 8.0, 99.0, 38.0, 43.0, 97.0, 82.0, 98.0, 60.0, 85.0, 75.0, 56.0, 70.0, 98.0, 64.0, 97.0, 43.0, 78.0, 81.0, 76.0, 7.0, 69.0, 100.0, 56.0, 71.0, 78.0, 66.0, 98.0, 86.0, 77.0, 39.0, 93.0, 51.0, 71.0, 100.0, 75.0, 57.0, 88.0, 61.0, 98.0, 65.0, 12.0, 4.0, 7.0, 41.0, 53.0, 77.0, 77.0, 75.0, 75.0, 62.0, 26.0, 94.0, 70.0, 94.0, 70.0, 58.0, 71.0, 17.0, 49.0, 71.0, 70.0, 60.0, 14.0, 55.0, 54.0, 32.0, 57.0, 80.0, 62.0, 79.0, 58.0, 81.0, 66.0, 69.0, 59.0, 66.0, 79.0, 46.0, 42.0, 67.0, 53.0, 47.0, 68.0, 68.0, 58.0, 47.0, 67.0, 51.0, 26.0, 34.0, 51.0, 57.0, 61.0, 62.0, 26.0, 37.0, 52.0, 48.0, 47.0, 47.0, 36.0, 53.0, 5.0, 38.0, 7.0, 22.0, 8.0, 26.0, 26.0, 20.0, 13.0, 41.0, 38.0, 49.0, 42.0, 39.0, 43.0, 40.0, 20.0, 33.0, 33.0, 33.0, 42.0, 32.0, 36.0, 39.0, 37.0, 36.0, 15.0, 37.0, 37.0, 22.0, 32.0, 21.0, 3.0, 32.0, 15.0, 22.0, 30.0, 27.0, 20.0, 23.0, 15.0, 19.0, 24.0, 20.0, 22.0, 4.0, 1.0, 11.0, 18.0, 21.0, 16.0, 10.0, 18.0, 11.0, 12.0, 13.0, 11.0, 15.0, 8.0, 16.0, 12.0, 11.0, 12.0, 15.0, 15.0, 6.0, 9.0, 9.0, 9.0, 6.0, 6.0, 5.0, 1.0, 6.0, 5.0, 4.0, 2.0, 2.0, 2.0], [1.0, 1.0, 1.0, 2.0, 1.0, 5.0, 3.0, 5.0, 5.0, 6.0, 4.0, 8.0, 8.0, 3.0, 9.0, 9.0, 5.0, 9.0, 10.0, 6.0, 8.0, 10.0, 8.0, 15.0, 12.0, 12.0, 4.0, 10.0, 17.0, 3.0, 16.0, 21.0, 17.0, 19.0, 27.0, 19.0, 23.0, 20.0, 9.0, 29.0, 14.0, 19.0, 29.0, 11.0, 6.0, 15.0, 23.0, 32.0, 7.0, 31.0, 23.0, 12.0, 34.0, 39.0, 32.0, 32.0, 12.0, 36.0, 19.0, 35.0, 28.0, 21.0, 16.0, 43.0, 27.0, 8.0, 17.0, 22.0, 34.0, 38.0, 39.0, 40.0, 42.0, 38.0, 46.0, 12.0, 30.0, 46.0, 45.0, 49.0, 11.0, 31.0, 48.0, 48.0, 54.0, 41.0, 45.0, 44.0, 49.0, 48.0, 31.0, 56.0, 19.0, 57.0, 48.0, 16.0, 33.0, 53.0, 43.0, 33.0, 40.0, 42.0, 62.0, 44.0, 57.0, 25.0, 4.0, 21.0, 59.0, 60.0, 68.0, 39.0, 59.0, 69.0, 30.0, 11.0, 25.0, 45.0, 45.0, 6.0, 75.0, 29.0, 33.0, 76.0, 65.0, 78.0, 48.0, 68.0, 60.0, 45.0, 57.0, 81.0, 53.0, 81.0, 36.0, 66.0, 69.0, 65.0, 6.0, 60.0, 87.0, 49.0, 63.0, 70.0, 60.0, 90.0, 79.0, 71.0, 36.0, 86.0, 48.0, 68.0, 96.0, 72.0, 55.0, 85.0, 59.0, 96.0, 64.0, 12.0, 4.0, 7.0, 41.0, 53.0, 78.0, 78.0, 76.0, 76.0, 63.0, 27.0, 98.0, 73.0, 99.0, 74.0, 62.0, 78.0, 19.0, 55.0, 80.0, 79.0, 68.0, 16.0, 63.0, 62.0, 37.0, 66.0, 95.0, 74.0, 95.0, 70.0, 98.0, 80.0, 84.0, 72.0, 81.0, 100.0, 60.0, 55.0, 88.0, 70.0, 64.0, 93.0, 95.0, 82.0, 67.0, 96.0, 74.0, 38.0, 50.0, 76.0, 85.0, 93.0, 95.0, 41.0, 59.0, 84.0, 81.0, 82.0, 83.0, 64.0, 95.0, 9.0, 70.0, 13.0, 41.0, 15.0, 49.0, 50.0, 39.0, 26.0, 83.0, 77.0, 100.0, 86.0, 83.0, 92.0, 87.0, 44.0, 73.0, 77.0, 78.0, 100.0, 78.0, 89.0, 98.0, 94.0, 92.0, 39.0, 99.0, 100.0, 62.0, 92.0, 61.0, 9.0, 98.0, 47.0, 70.0, 99.0, 93.0, 70.0, 81.0, 55.0, 72.0, 93.0, 78.0, 86.0, 16.0, 4.0, 45.0, 77.0, 92.0, 73.0, 46.0, 86.0, 54.0, 60.0, 65.0, 55.0, 80.0, 43.0, 86.0, 67.0, 65.0, 76.0, 98.0, 99.0, 44.0, 73.0, 73.0, 77.0, 69.0, 76.0, 68.0, 14.0, 84.0, 73.0, 87.0, 71.0, 71.0, 76.0], [28, 27, 27, 27, 26, 24, 24, 22, 22, 22, 21, 18, 17, 17, 16, 16, 14, 13, 13, 10])

instances = [instance1,instance2,instance3,instance4,instance5,instance6,instance7,instance8,instance9,instance10,instance11,instance12]

for test in instances
     println(@elapsed branchandbound_HungFisk( test.c, test.w, test.W)[2])
end

