#= Optimisation Discrète et Combinatoire
M01KP problem =#

using JuMP, GLPK

function modelM01KP(n::Int64, m::Int64, couts::Vector{Int64}, poids::Vector{Int64}, capa::Vector{Int64})

    model::Model = Model(GLPK.Optimizer)

    @variable(model, x[1:m, 1:n] >= 0)

    @objective(model, Max, sum(sum(couts[j] * x[i, j] for j = 1:n) for i = 1:m))

    @constraint(model, sum(poids[j] * x[1, j] for j = 1:n) <= capa[1])
    @constraint(model, sum(poids[j] * x[2, j] for j = 1:n) <= capa[2])
    @constraint(model, cst[j = 1:n], sum(x[i, j] for i = 1:m) <= 1)

    optimize!(model)

    # Affichages
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x)
end

function solve_modelM01KP_surrogate(nb_objets::Int64, nb_sacs::Int64, couts::Vector{Float64}, poids::Vector{Float64}, capa::Vector{Int64}, coeff::Float64, S::Vector{Int64})

    model::Model = Model(GLPK.Optimizer)

    @variable(model, x[1:(nb_sacs + 1), 1:nb_objets], Bin)

    @objective(model, Max, sum(sum(couts[j] * x[i, j] for j = 1:nb_objets) for i = 1:nb_sacs))

    @constraint(model, sum(coeff * sum(poids[j] * x[i, j] for j = 1:nb_objets) for i = 1:nb_sacs) <= sum(coeff * capa[i] for i = 1:nb_sacs))
    @constraint(model, cst[j = 1:nb_objets], sum(x[i, j] for i = 1:(nb_sacs + 1)) <= 1)
    @constraint(model, cstbis[j in S], sum(x[i, j] for i = 1:(nb_sacs + 1), j in S) == 1)
    optimize!(model)

    return value.(x), objective_value(model)
end

function solveM01KP()
    n::Int64 = 6
    m::Int64 = 2
    couts::Vector{Int64} = [110, 150, 70, 80, 30, 5]
    poids::Vector{Int64} = [40, 60, 30, 40, 20, 5]
    capa::Vector{Int64} = [65, 85]

    modelM01KP(n, m, couts, poids, capa)
end
#=
function solveM01KP_surrogate()
    n::Int64 = 6
    m::Int64 = 2
    couts::Vector{Int64} = [110, 150, 70, 80, 30, 5]
    poids::Vector{Int64} = [40, 60, 30, 40, 20, 5]
    capa::Vector{Int64} = [65, 85]

    modelM01KP_surrogate(n, m, couts, poids, capa, 2)
end
=#
