__precompile__(true)

module CorrelatedEquilibria
#= Functions
    build_game(path): builds the game based on the struct with the data from yaml file
    maxPayoff(game): finds the maximum Payoff that exists for any player in a game
    minPayoff(game): finds the minimum Payoff that exists for any player in a game
    initializeV(g, gamma): initializes the Hypercube at the start of the algorithm
    computeQ(g, V, gamma): Compute all value vectors 
    updatePayoffs(Q, g, gamma, outerError=0.0001): performs one step of the algorithm
    computeCE(g, gamma, steps=3, plotting=True): this is the function with executes the algorithm as often as the steps indicate
=#

include("Utils.jl")
include("Plotting.jl")
using LazySets
using Polyhedra
using LinearAlgebra
using Plots
using CDDLib
using JuMP

# construct the struct for the MarkovGame
struct MarkovGame
    name::String
    players::Vector{Any}
    states::Vector{Any}
    actions::Dict{Any, Any}
    joint_actions::Dict{Any, Any}
    payoffs::Dict{Any, Any}
    probabilities::Dict{Any, Any}
end

function build_game(path)
    #= builds the game based on the struct with the data from yaml file =#
   
    # read the file
    data = read_yaml(path)

    # checking the yaml file for the correct structure
    check_yaml(data)
    
    game = MarkovGame(
        data["name"],
        data["players"],
        data["states"],
        data["actions"],
        data["joint_actions"],
        data["payoffs"],
        data["probabilities"])

    return game
end

function maxPayoff(game) 
    #= finds the maximum Payoff that exists for any player in a game
    Params:
        game: Markov Game
        idx: Array with the indices of all joint actions
    =#
    return extremal(game, "max")
end

function minPayoff(game) 
    #= finds the minimum Payoff that exists for any player in a game
    Params:
        game: Markov Game
        idx: Array with the indices of all joint actions
    =#
    return extremal(game, "min")
end

function initializeV(g, gamma)
    #= initializes the Hypercube at the start of the algorithm
    Params:
        g (Struct Markov Game): Stochastic Game
        gamma (Float): discount factor 
    Returns:
        V (Vector): Inital Hypercube
    =#
    dimension = length(g.players)
    In = Matrix{Float64}(1.0*I, dimension, dimension)
    min = minPayoff(g)/(1-gamma)
    max = maxPayoff(g)/(1-gamma)

    # initialise Hypercube
    V::AbstractArray{LazySet,1} = [
          HPolyhedron([
                  LinearConstraint(k*In[:,i] ,k*extr) 
                  for i in 1:dimension for (k,extr) in [(-1,min),(1,max)]       
                          ]) 
            for s in g.states]
end

function computeQ(g, V, gamma)
    #= Compute all value vectors 
    Params:
        g (Struct Markov Game): Stochastic Game
        V (Vector): Inital Hypercube
        gamma (Float): discount factor 
    Returns:
        Q (Vector): Contains all value vectors 
    =#
    model = Model(CDDLib.Optimizer{Float64})

    # directPayoffs = gamma * sum(P(nextState | (state, jointAction) * V(nextState)))
    # thearray = R(state, jointAction)
    Q = [   [
            let directPayoffs = [g.payoffs[s][jointActID][pID]
                    for (pID, p) in enumerate(g.players)]
                thearray = 
                vcat([Singleton(directPayoffs)],[
                        gamma * g.probabilities[s][jointActID][nextS] * V[nextID]
                        for (nextID, nextS) in enumerate(g.states)])
                MinkowskiSumArray(thearray)
                end
            for jointActID in 1:length(g.payoffs[s])
            ]   
        for s in g.states
        ]
    return Q
end

function updatePayoffs(Q, g, gamma, outerError)
    #= performs one step of the algorithm
    Params:
        Q (Vector): Contains all value vectors 
        g (Struct MarkovGame): Stochastic game
        gamma (Float): discount factor
        outerError (Float): outer Error for over approximation
    Returns:
        V (Vector): Inital Hypercube
        X (Matrix): Constraint matrix
    =#
    minima = computeQsMinima(Q)
    actionsArray = actionsArrays(g)

    N = length(g.players)
    V = Array{LazySet,1}(undef,length(g.states))
      
    stateCount = length(Q)
    if stateCount ≠ length(g.states)
      print("Sth. wrong with the number of states")
    end  
    X = Array{LazySet,1}(undef,stateCount)

    for (sID,s) in enumerate(g.states)
        actionCount = length(g.joint_actions[s])

        # the total dimension of the new polytope will be 
        total = N + actionCount*N + actionCount
        It = Matrix{Float64}(1.0*I,total,total)

        function equation19()
            # EQUATION 19: v - sum(qbar[state][joint action]) = 0
            # Note: equation of inequality: x1-x2=0 if x1-x2 <= 0 and -x1+x2 <= 0

        
            X_tmp = let theConstr = [
                LinearConstraint( k*vOf(It, pID) - sum(k * [qbarOf(It, N, pID, jointActID) for jointActID in 1:length(keys(g.payoffs[s]))]), 0.) 
                for (pID,p) in enumerate(g.players) for k in [1,-1]
                    ]
            LazySets.HPolytope(theConstr)
            end
            return X_tmp
        end
        

        function equation20(s, minima)
            # 0 >= (sum(omega[joint action] * Qbar[state][deviated joint action] for all possible deviations) - Q[state][joint action]

            for (pID, p) in enumerate(g.players)
        
                for proposed in 1:length(g.actions[s][p])
                    # save joint action with proposed action for player p
                    compatible = [ jointActTup for jointActTup in actionsArray[sID] if jointActTup[pID]==proposed ]
                    # compute the sum of matrix rows
                    lhs = sum([qbarOf(It, N, pID, findfirst(isequal(jointActTup), actionsArray[sID])) for jointActTup in compatible])
                    # check for different actions
                    for other in 1:length(g.actions[s][p])
                        if other == proposed
                            continue
                        end

                        # compute maximum possible punishment
                        rhs = sum([ minima[sID][findfirst(isequal(subst(jointActionTup, pID, other)), actionsArray[sID])][pID] * omegaOf(It, N, actionCount, findfirst(isequal(jointActionTup), actionsArray[sID]))
                        for jointActionTup in compatible
                        ])
                    
                        addconstraint!(X[sID], LinearConstraint(rhs-lhs,0.))
                    end
                end
            end
        end

        function equation21()
            # newvector + scalar * omega[joint action] 
            for jointActID in 1:actionCount
            
                pre = jointActID*N 
                post = total - (pre + N)
                oldQas = overapproximate(Q[sID][jointActID], HPolygon, outerError)
                oldConstraints = constraints_list(oldQas)
                for c in oldConstraints
                    vector::Array{Float64,1} = c.a
                    scalar::Float64 = c.b
                    prefix = zeros(Float64, pre)
                    postfix = zeros(Float64, post)
                    newvector = vcat(prefix, vector, postfix) - scalar*omegaOf(It, N, actionCount, jointActID)
                    
                    addconstraint!(X[sID], LinearConstraint(newvector,0.))
                end
            end 
        end

        function equation22()
            # sum(omega) = 1
            for k in [1.0, -1.0] 
                addconstraint!(X[sID], LinearConstraint( sum([k*omegaOf(It, N, actionCount, act) for act in 1:length(keys(g.payoffs[s]))]), k))
              end 
        end

        function equation23()
            # omega(joint action) >= 0
            for act in 1:length(keys(g.payoffs[s]))
                addconstraint!(X[sID], LinearConstraint(-omegaOf(It, N, actionCount, act),0.))
              end
        end
        
        ### EQUATION 19:
        print("now equation ", 19, "...") 
        X[sID] = equation19()
        println("done.")

        ### EQUATION 2O
        print("now equation ", 20, "...")
        equation20(s, minima)
        println("done.")

        ### EQUATION 21
        print("now equation ", 21, "...")
        equation21()
        println("done.")
        
        ### EQUATION 22
        print("now equation ", 22, "...")
        equation22()
        println("done.")

        ### EQUATION 23
        print("now equation ", 23, "...")
        equation23()
        println("done.")

        proj = Matrix{Float64}(1.0*I,N,total)
        V[sID] = LazySets.LinearMap(proj,X[sID])
    end

    return (V,X)
end

function computeCE(g, gamma, steps=3, plotting=true, outerError=0.0001)
    #= this is the function with executes the algorithm as often as the steps indicate
    Params:
        g (Struct Markov Game): Stochastic Game
        gamma (Float): discount factor 
        steps (Int): number of steps the algorithm should be executed
        plotting (Bool): if gif should be produced or not
        outerError (Float): outer Error for over approximation
    Returns:
        V (Vector): Set of value vectors that are in the correlated equilibria
        X (Matrix): Constraint vectors for final V
        animArr (gif): Set of plots that show progress of V for each state
    =#
    # init
    V = initializeV(g,gamma)
    oldminis = nothing
    minis = nothing

    animArr = Array{Any}(undef, length(g.states))
    for (sID, s) in enumerate(g.states)
        animArr[sID] = Plots.Animation()
    end

    # algorithm loop   
    for k in 1:steps
        # append plot
        if plotting
            for (sID, s) in enumerate(g.states)
                plot(V[sID], title="V[$sID]: Step $k")
                Plots.frame(animArr[sID])
            end
        end
        
        # perform one step of algorithm
        println("update step ", k)
        Q = computeQ(g, V, gamma)
        oldminis = minis        
        minis = computeQsMinima(Q)
        (V,X) = updatePayoffs(Q, g, gamma, outerError)
      
        # termination condition
        # TODO: termination depending on update size
        if k == steps
            return (V,X, animArr)
        end
    end
end

# if you export a function the pkg name does not need to be written in front of it when executing the function
# don't export internal functions (only ones for public use)
export build_game
export maxPayoff
export minPayoff
export computeQ
export initializeV
export computeCE

export plot_V_progress
export plot_V

end # end of module