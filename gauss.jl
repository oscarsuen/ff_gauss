using Statistics
using StatsPlots

# Q_DEFAULT = 2
# DIMS_DEFAULT = (5, 3)
# SAMPLES_DEFAULT = 100

# TODO: type + prime checking
struct FF{q}
    x::Int
    FF{q}(x::Int) where {q} = new{q}(mod(x, q))
end
Base.show(io::IO, a::FF{q}) where {q} = print(io, a.x, "q", q)
Base.Broadcast.broadcastable(a::FF) = Ref(a)

Base.zero(::Type{FF{q}}) where {q} = FF{q}(0)
Base.zero(::FF{q}) where {q} = FF{q}(0)
Base.one(::Type{FF{q}}) where {q} = FF{q}(1)
Base.one(::FF{q}) where {q} = FF{q}(1)

function make_qinv(q::Int)
    qinv = zeros(Int, q-1)
    for i ∈ 1:(q-1)
        for j ∈ 1:(q-1)
            if mod(i*j, q) == 1
                qinv[i] = j
                continue
            end
        end
    end
    # println("made qinv $q")
    # display(qinv)
    return qinv
end
qinv_dict = Dict{Int, Vector{Int}}()

function inv(a::FF{q}) where {q}
    qinv = get!(qinv_dict, q) do 
        make_qinv(q)
    end
    if iszero(a)
        # throw(DivideError(x, "cannot invert 0"))
        error("cannot invert 0")
    end
    return FF{q}(qinv[a.x])
end

Base.:+(a::FF{q}, b::FF{q}) where {q} = FF{q}(a.x+b.x)
Base.:*(a::FF{q}, b::FF{q}) where {q} = FF{q}(a.x*b.x)
Base.:-(a::FF{q}) where {q} = FF{q}(-a.x)
Base.:-(a::FF{q}, b::FF{q}) where {q} = a + (-b)
Base.:/(a::FF{q}, b::FF{q}) where {q} = a * inv(b)

gen_matrix(q::Int, dims::Tuple{Int, Int}) = FF{q}.(rand(0:(q-1), dims...))

function gelim!(A::Matrix)
    # Does not back substitute
    G = A
    # G = copy(A)
    # display(G)
    i, j = 1, 1
    r, c = size(G)
    while i ≤ r && j ≤ c
        f = findfirst(!iszero, G[i:end, j])
        if isnothing(f)
            j += 1
            continue
        end

        f = i + f-1
        row_i = copy(G[i, :])
        row_f = copy(G[f, :])
        G[i, :] = row_f
        G[f, :] = row_i
        # println("swapped rows $i $f")
        # display(G)

        @assert !iszero(G[i, j])
        G[i, :] ./= G[i, j]
        # println("normalized row $i")
        # display(G)

        for k ∈ (i+1):r
            if !iszero(G[k, j])
                G[k, :] -= G[k, j] .* G[i, :]
                # println("subtracted $(G[k, j]) times row $i from row $k")
                # display(G)
            end
        end

        i += 1
        j += 1
    end
    return G
end

function isrowechelon(G::Matrix)
    # println("check matrix")
    # display(G)
    j = 0
    r, c = size(G)
    for i ∈ 1:r
        f = findfirst(!iszero, G[i, :])
        # display("i=$i, j=$j, f=$f")
        if isnothing(j)
            if !isnothing(f)
                # println("row $i should be all 0s")
                return false
            end
        elseif !isnothing(f)
            if !isone(G[i, f])
                # println("row $i doesn't have leading 1")
                return false
            end
            if f ≤ j
                # println("row $i misshaped")
                return false
            end
        end
        j = f
    end
    return true
end

# function rankgauss(G::Matrix)
#     return count(!iszero, eachrow(G))
#     return sum((!iszero).(eachrow(G)))
#     return searchsortedfirst(1:size(G)[1], 1, by=x->iszero(G[x, :])) #fails
#     return searchsortedfirst(eachrow(G), 0, by=iszero)-1
# end
rankgauss(G::Matrix) = searchsortedfirst(eachrow(G), 0, by=iszero) - 1

# function rank!(A::Matrix)
#     display(A)
#     G = gelim!(A)
#     display(G)
#     return rankgauss(G)
# end
rank!(A::Matrix) = rankgauss(gelim!(A))

rank_distribution(q::Int, dims::Tuple{Int, Int}, samples::Int) = (rank!(gen_matrix(q, dims)) for _ ∈ 1:samples)

square_ranks(q::Int, samples::Int, max_k::Int) = (rank_distribution(q, (k, k), samples) for k ∈ 1:max_k)
mean_ranks(q::Int, samples::Int, max_k::Int) = mean.(square_ranks(q, samples, max_k))

function gen_linegraph(qs::AbstractVector{Int}, ns::AbstractVector{Int}, samples::Int; diff=true::Bool)
    if diff
        plot()
    else
        plot(ns, ns, label="0")
    end
    for q ∈ qs
        if diff
            errorline!(ns, stack(n .- rank_distribution(q, (n, n), samples) for n ∈ ns), errorstyle=:stick, errortype=:sem, secondarycolor=:matched, label="$q")
        else
            errorline!(ns, stack(rank_distribution(q, (n, n), samples)|>collect for n ∈ ns), errorstyle=:stick, errortype=:sem, secondarycolor=:matched, label="$q")
        end
    end
    # errorline(ns, stack(stack(n .- rank_distribution(q, (n, n), samples) for n ∈ ns; dims=1) for q ∈ qs), errorstyle=:stick, errortype=:sem, label=qs')
    gui()
end

# histogram(rank_distribution(2, (100, 100), 1000)|>collect); gui()
