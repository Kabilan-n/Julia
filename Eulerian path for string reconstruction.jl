A = "GATTTTCCCACGGCCACGCCCGGATATGAGGTAATTTTGGGCGGTTGCAAATAAAATTAGGACATGGTGGCG"
function spliting(A,k)
    X = length(A)-k+1
    B =[]
    for i in 1:X
        D = A[i:i+k-1]
        append!(B,[D])
    end
    return B
end
B=spliting(A,3)
print(B)
length(B)

function del!(a, item)
    q = findall(x->x==item, a)
    deleteat!(a,q[1])
end

function sorting(B)
    k = length(B[1])
    R = Dict("A" => 1, "C" =>2, "G" => 3, "T" =>4)
    O = Dict()
    no =[]
    Sorted = []
    for j in 1:length(B)
        Y = B[j]
        n=0
        for i in 1:length(Y)
            x  =R[string(Y[i])]*(10^(k-i))
            n = x + n
        end
        O[n]=Y
        append!(no,n)
    end
    for i in 1:length(no)
        a =minimum(no)
        del!(no,a)
        append!(Sorted,[O[a]])
    end
    
    
    return Sorted
end
B=sorting(B)
print(B)
length(B)

function count(list)
    frequency = Dict()
    for item in list
        if haskey(frequency, item)
            frequency[item] += 1
        else
            frequency[item] = 1
        end
    end
    return frequency
end


using DataStructures
function BruijnKmers(text)
    prefix =[]
    suffix =[]
    dictionary = Dict()
    indegree = Dict()
    outdegree = Dict()
    for element in text
        kmer = element[1:length(element)-1]
        kmer2 = element[2:length(element)]
        append!(prefix,[kmer])
        append!(suffix,[kmer2])
        if length(kmer2) == length(prefix[1])
            if haskey(dictionary, kmer)
                append!(dictionary[kmer],[kmer2])
            else
                dictionary[kmer] = [kmer2]
            end
        end
    end
    outd = count(prefix)
    ind = count(suffix)
    result = dictionary
    return result,ind,outd
end

result,ind,outd = BruijnKmers(B)
println(result)
println(ind)
println(outd)

function stacking(result,start)
    res = copy(result)
    s = Stack{String}() 
    circuit = Stack{String}()
    push!(s, start)
    current = start
    a =1
    while length(s)>0
        if current in keys(result) && length(result[current]) > 0
            next =res[current]
            next = next[length(next)]
            push!(s, next)
            X = res[current]
            del!(X,next)
            res[current] = X
            current = next

        else
            w = pop!(s)
            push!(circuit,w)
            current = w
        end
    end
    return circuit
end
circuit = stacking(result,"GA")

#println(circuit)



function seq(que)
    x = Int(length(first(que)))
    z = length(que)
    final = string()
    for i in 1:length(que)
        q = pop!(que)
        final = string(final,string(q[1]))
        if i == z
            final = string(final,string(q[2:x]))
        end
    end
    return final
end
a = seq(circuit)



"GATTTTCCCACGGCCACGCCCGGATATGAGGTAATTTTGGGCGGTTGCAAATAAAATTAGGACATGGTGGCG"
