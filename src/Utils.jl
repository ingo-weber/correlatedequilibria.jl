import YAML

#= Functions:
    read_yaml(path): reads a yaml file (.yml or .yaml) into a data variable
    check_yaml(data): builds a game based on the defined struct using the data from the yaml file
    extremal(g, type): calculating the extremal payoffs of a Markov game (auxiliary function) 
    computeMinima(P): computes the minimum for each player in present state
    computeQsMinima(Q): return the minimum of each Q[state][jointAction]
    actionsArrays(g): return an array with the indvidual action indices of each player in every joint action 
    vOf(It, pID): returns the vector of It of the given dimension
    qbarOf(It, N, pID, jointActID): returns the vector of It for a specific joint action and player
    omegaOf(It, N, actionCount, jointActID): returns the vector of It for the omega of a specific joint action
    indexOf(g, s, pID, a): returns the index of an action of a player
    subst(jointActionTup, pID, other): substitutes the deviated action of a player for proposed joint action
=#

function read_yaml(path)
    #= reads a yaml file (.yml or .yaml) into a data variable =#
    if occursin(".yml", path) || occursin(".yaml", path)
        try 
            data = YAML.load_file(path)
            return data
        catch e
            throw(e)
        end
    else
        println("Wrong file format; file name needs to end with .yml or .yaml")
    end
end

function check_yaml(data)
    #= builds a game based on the defined struct using the data from the yaml file =#
    
    # checking the name
    try
        data["name"]
    catch e
        println("The yaml file is missing the name of the game")
        throw(e)
    end

    if typeof(data["name"]) != String
        println("Please enter a String for the file name")
    end

    # checking the players
    try
        data["players"]
    catch e
        println("Please enter the players of the game")
        throw(e)
    end

    #=
    if typeof(data["players"]) != Vector{Any}
        println("Wrong format! Pleaser enter the players as a list of Strings like: ")
        println("Players:\n   - Player1\n   - Player2")
    end
    =#

    # TODO: checking the rest of the file
    # TODO: check that all joint action keys are consecutive ints (1,2,3,4,5..)
end

function extremal(g, type)
#= calculating the extremal payoffs of a Markov game (auxiliary function) 
    Params:
        g (struct MarkovGame): Stochastic game
        type (String): type of desired extreme (max or min)
    Returns:
        extreme (Float): extremum
    =#

    # init return variables
    extreme = nothing

    # iterate through state 
    for s in keys(g.payoffs)
        # iterate throug joint actions
        for ja in keys(g.payoffs[s])
            # iterate to individual rewards
            for rew in g.payoffs[s][ja]
                # if extreme has not been initialised
                if extreme === nothing
                    extreme = rew
                # if the maximum is required to be found
                elseif type == "max"
                    if extreme < rew
                        extreme = rew
                    end
                # if the minimum is required to be found
                elseif type == "min"
                    if extreme > rew
                        extreme = rew
                    end
                end
            end
        end
    end
    return extreme
end

function computeMinima(P)
    #= computes the minimum for each player in present state
    Params:
        P (LazySets.Polyhedra): Q[state][jointAction]
    Returns:
        the support function of a set in a given direction
        
    =#
    n = LazySets.dim(P)
    In = Matrix(I, n, n)
    concreteP = P
     return [-LazySets.ρ(Vector(-In[:,i]),concreteP)  for i in 1:n]
end

function computeQsMinima(Q)
    #= return the minimum of each Q[state][jointAction]
    Params:
        Q (Vector): Contains all value vectors 
    Returns:
        Qs (Vector): Vector with all the minimas
    =#
    Qs = [ 
      [
       computeMinima(Q[s][coord])
        for coord in CartesianIndices(Q[s])
        ]
      for s in 1:length(Q)
    ]
    return Qs 
end

function actionsArrays(g)
    #= return an array with the indvidual action indices of each player in every joint action 
    Params:
        g (Struct MarkovGame): Stochastic game
    Returns:
        actsArr (Array): holds CartesianIndices of each joint action
    =#

    # construct empty array
    actsArr = [ Array{Any}(undef, length(g.joint_actions[s])) for s in g.states]

    for (sID, s) in enumerate(g.states)
        
        for (jointActKey, jointAct) in g.joint_actions[s]
            # build a tuple with the indices of the individual action of a player
            # for each joint action e.g.: 
            # jointAct1 = ["Action1" "Action2"] --> (1 2)
            tup = Tuple(indexOf(g, s, pID, a) for (pID, a) in enumerate(jointAct))
            actsArr[sID][jointActKey] = CartesianIndex(tup)
        end
    end
    return actsArr
end

function vOf(It, pID)
    #= returns the vector of It of the given dimension
    Params:
        It (Matrix): holds indices for constraints
        pID (int): number of player
    =#
    It[:, pID]
end

function qbarOf(It, N, pID, jointActID)
    #= returns the vector of It for a specific joint action and player
    Params:
        It (Matrix): holds indices for constraints
        N (int): total number of players
        pID (int): number of player
        jointActID (int): number of joint action
    =#
    It[:, ((jointActID*N) + pID)]
end

function omegaOf(It, N, actionCount, jointActID)
    #= returns the vector of It for the omega of a specific joint action
    Params:
        It (Matrix): holds indices for constraints
        N (int): total number of players
        actionCount (int): total number of actions
        jointActID (int): number of joint action
    =#
    It[:, (N + N*actionCount + jointActID)]
end

function indexOf(g, s, pID, a) 
    #= returns the index of an action of a player
    Params:
        g (Struct MarkovGame): Stochastic game
        s (String): state name
        pID (Int): player number
        a (String): name of individual action
    =#
    return findfirst(isequal(a), g.actions[s][g.players[pID]])
end

function subst(jointActionTup, pID, other)
    #= substitutes the deviated action of a player for proposed joint action
    Params: 
        jointActionTup (CartesianIndex): tuple with indices of individual player actions of joint action
        pID (Int): player number of player who deviates
        other (Int): index of deviated player actionCount
    Returns:
        CartesianIndex of deviated joint action
    =#
    a = [x for x in Tuple(jointActionTup)]
    a[pID] = other
    CartesianIndex(Tuple(a))
end