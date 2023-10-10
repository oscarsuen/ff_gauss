using Statistics

q = 2
n = 3
m = 3

struct FF
    x::Int
    FF(x) = new(mod(x, q))
end
Base.show(io::IO, a::FF) = print(io, a.x, "q", q)
Base.Broadcast.broadcastable(a::FF) = Ref(a)

Base.zero(::Type{FF}) = FF(0)
Base.zero(::FF) = FF(0)
Base.one(::Type{FF}) = FF(1)
Base.one(::FF) = FF(1)

qinv = zeros(Int, q-1)
for i ∈ 1:(q-1)
    for j ∈ 1:(q-1)
        if mod(i*j, q) == 1
            qinv[i] = j
            continue
        end
    end
end
function inv(a::FF)
    if iszero(a)
        # throw(DivideError(x, "cannot invert 0"))
        error("cannot invert 0")
    end
    return FF(qinv[a.x])
end

Base.:+(a::FF, b::FF) = FF(a.x+b.x)
Base.:*(a::FF, b::FF) = FF(a.x*b.x)
Base.:-(a::FF) = FF(-a.x)
Base.:-(a::FF, b::FF) = a + (-b)
Base.:/(a::FF, b::FF) = a * inv(b)

gen_matrix(rows=n, cols=m) = FF.(rand(0:(q-1), rows, cols))

function gelim!(A::Matrix{FF})
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

function isrowechelon(G::Matrix{FF})
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

function rankgauss(G::Matrix)
    # return count(!iszero, eachrow(G))
    # return sum((!iszero).(eachrow(G)))
    # return searchsortedfirst(1:size(G)[1], 1, by=x->iszero(G[x, :])) #fails
    return searchsortedfirst(eachrow(G), 0, by=iszero)-1
end
function rank!(A::Matrix{FF})
    return rankgauss(gelim!(A))
    # display(A)
    # G = gelim!(A)
    # display(G)
    # return rankgauss(G)
end

function rank_distribution(samples, rows=n, cols=m)
    return (rank!(gen_matrix(rows, cols)) for _ in 1:samples)
end

