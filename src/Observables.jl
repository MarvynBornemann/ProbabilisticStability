function get_convergence_to_state(sol::AbstractODESolution, state, distance; tail_frac=0.8)
    L = length(sol.t)
    idx = max(1, round(Int, tail_frac*L)):L
    x = sol.t[idx]
    y = [distance(state, p) for p in sol.u[idx]]
    X = zeros(length(x),2)
    X[:,1] = x
    X[:,2] .= 1.0
    slope, intercept = X \ log.(y .+ 1E-10) # (y .+ 1E-10) to avoid log(0.0) if y is exactly 0.0 somewhere
    return slope
end

function eval_convergence_to_state(sol::AbstractODESolution, state, distance; tail_frac=0.8, verbose=false)
    slope = get_convergence_to_state(sol, state, distance; tail_frac=tail_frac)
    if verbose
        println("The estimated slope is $slope.")
    end
    return slope < 0
end

function get_final_distance_to_state(sol::AbstractODESolution, state, distance; threshold=1E-3, verbose=true)
    return distance(state, sol[end])
end

function eval_final_distance_to_state(sol::AbstractODESolution, state, distance; threshold=1E-3, verbose=true)
    d = distance(state, sol[end])
    if verbose
        println("The final state distance is $d.")
    end
    return d < threshold
end

function eval_final_distance_to_state(d::Real; threshold=1E-3, verbose=true)
    if verbose
        println("The final state distance is $d.")
    end
    return d < threshold
end

"""
e.g. distance -> PeriodicEuclidean([Inf,2π]) or Euclidean()
"""
function eval_final_distance_to_state(pint, state, distance; threshold=1E-3, verbose=true)
    d = colwise(distance, state, pint.u)
    if verbose
        println("The final state distance is in the range: ", extrema(d))
    end
    return d .< threshold
end

# function eval_final_distances_to_state(sol, state, distance; state_filter=nothing, threshold=1E-3, verbose=true)
#     d = distance.(state, sol[end])
#     if isnothing(state_filter)
#         dmax = maximum(d)
#     else
#         dmax = maximum(d[state_filter])
#     end
#     if verbose
#         println("The max final state distance across all dimensions is $dmax.")
#     end
#     return dmax < threshold
# end

function get_mean_distance_to_state(sol::AbstractODESolution, state, distance; tail_frac=0.8)
    L = length(sol.t)
    idx = max(1, round(Int, tail_frac*L)):L
    return mean([distance(state, p) for p in sol.u[idx]])
end

function eval_mean_distance_to_state(sol::AbstractODESolution, state, distance; threshold=1E-3, tail_frac=0.8, verbose=true)
    d = get_mean_distance_to_state(sol, state, distance; tail_frac=tail_frac)
    if verbose
        println("The mean state distance is $d.")
    end
    return d < threshold
end

eval_mean_distance_to_state(d; threshold=1E-3) = d < threshold

function get_max_distances_to_state(sol::AbstractODESolution, state, distance; tail_frac=0.8)
    L = length(sol.t)
    #idx = max(1, round(Int, tail_frac*L)):L
    return maximum(distance.(fp, sol), dims=2)
end

function eval_max_distances_to_state(sol::AbstractODESolution, state, lb, ub, distance; tail_frac=0.8, verbose=true)
    d = get_max_distances_to_state(sol, state, distance; tail_frac=tail_frac)
    if verbose
        println("The max state distances are $d.")
    end
    return lb .< d .< ub
end

eval_max_distances_to_state(d, lb, ub) =  lb .< d .< ub

function eval_max_distance_to_state(sol::AbstractODESolution, state, lb, ub, distance; tail_frac=0.8, verbose=true)
    d = get_max_distances_to_state(sol, state, distance; tail_frac=tail_frac)
    if verbose
        println("The max state distance is $d.")
    end
    return all(lb .< d .< ub)
end

eval_max_distance_to_state(d, lb, ub) =  all(lb .< d .< ub)

function get_trajectory_within_bounds(sol::AbstractODESolution, lb, ub; tail_frac=0, verbose=true)
    L = length(sol.t)
    idx = max(1, round(Int, tail_frac*L)):L
    return [all(lb .< p .< ub) for p in sol.u[idx]]
end

function eval_trajectory_within_bounds(sol::AbstractODESolution, lb, ub; tail_frac=0, verbose=true)
    # use Inf as bounds for dimension that should be beglected
    ep = get_trajectory_within_bounds(sol, lb, ub; tail_frac=tail_frac, verbose=verbose)
    #findlast(.! ep)
    return all(ep)
end

eval_trajectory_within_bounds(twb::Array{Bool,1}) = all(twb)

function all_evals(sol, i) # output_func
    a = eval_convergence_to_state(sol, fp, euclidean; verbose = false)
    b = eval_final_distance_to_state(sol, fp, euclidean; verbose = false)
    c = eval_trajectory_within_bounds(
        sol,
        [-2.0, -Inf],
        [2.0, Inf];
        verbose = false,
    )
    ([a; b; c], false)
end


"""
    limiting_curve(t_range)

Calculates a time series for a "low voltage fault-ride through limiting curve" based on the one given in:
"Probabilistic Stability Assessment for Active Distribution Grids" 
https://arxiv.org/abs/2106.09624
"""
function limiting_curve(t_range)
    slope = 0.7 / 2.85
    intercept = 0.15 * (1 - slope)

    limiting_curve = zeros(length(t_range))
    for i in 1:length(t_range)
        t = t_range[i]
        if t <= 0.15
            limiting_curve[i] = 0.15
        elseif t > 0.15 && t <= 3.0
            limiting_curve[i] = slope * t + intercept
        elseif t > 3.0 && t < 60.0
            limiting_curve[i] = 0.85
        elseif t >= 60.0
            limiting_curve[i] = 0.95
        end
    end
    return limiting_curve
end

function voltage_condition_surv(pg::PowerGrid, sol::AbstractODESolution)
    #catch case, where sol.t has less than 2 values
    if(length(sol.t) < 2)
        return false
    end

    limiting_curve_vol = limiting_curve(sol.t) # Generate limiting curve for this time series
    low_voltage_condition = Vector{Bool}(undef, length(pg.nodes))

    #my_sol = PowerGridSolution(sol, pg)
    my_sol = [State(pg, u) for u in sol.u]
    for n in 1:length(pg.nodes) # check all nodes
        v_node = []
        # for t in sol.t
        #     push!(v_node, my_sol(t, n, :v))
        # end
        for i in 1:length(sol.t)
            push!(v_node, my_sol[i][n, :v])
        end
        under_vol = findall(v_node .< limiting_curve_vol)
        if under_vol != Int64[]
            low_voltage_condition[n] = false # Condition is violated if a voltage is lower than in limiting curve
        else 
            low_voltage_condition[n] = true
        end
    end
    return all(low_voltage_condition) # Check if low voltage condition is not violated for all nodes
end

function eval_final_frequency(pg::PowerGrid, sol; threshold=0.18)
    #catch case, where sol.t has less than 2 values
    if(length(sol.t) < 2)
        return false
    end

    N = length(pg.nodes)
    f = zeros(N)
    t_last2 = sol.t[end-1:end]
    state_last2 = [State(pg, u) for u in sol.u[end-1:end]]
    for n in 1:N # check all nodes
        u_last2 = [state[n, :u] for state in state_last2] #last two complex voltages
        φ_diff = angle(u_last2[2] / u_last2[1]) #[rad] #u₂/u₁ = (v₂*exp(im*φ₂)) / (v₁*exp(im*φ₁)) = v₂/v₁ * exp(im*(φ₂-φ₁))
        t_diff = t_last2[2] - t_last2[1] #[s]
        ω = φ_diff / t_diff #[rad/s]
        f[n] = ω /(2π) #[Hz]
    end
    return all(abs.(f) .< threshold) #allowed change in frequency is +/- 180mHz
end