import Pkg
Pkg.add("GraphRecipes")


g = [0 1 1;
     0 0 1;
     0 1 0]

graphplot(g, names=1:3, curvature_scalar=0.1)

novertex = parse(Int,readline()) #getting no.of vertices
Dict1 = Dict()
Vertex =[]
E=[]
for i in 1:novertex
    append!(Vertex,readline())  #getting name of vertices
    Dict1[Vertex[i]]=i
end
print(Vertex)
print(Dict1)

#adding edges to the vertices
function Append(x)
    if (length( findall( y -> y == x, E )) == 0)
        append!(E,[x])
    end  
    
    return E
end

Append("ca") #ca
Append("cb") #cb
Append("cd") #cd
Append("ec") #ec
println(E)

# creates adjacency matrix
function Adjmat(Dict1,E)
    F = []
    T = []
    l = length(Dict1)
    A = zeros(Int64,(l,l))
    for j in 1:length(E)
        append!(F,E[j][1])
        append!(T,E[j][2])
        A[Dict1[F[j]] , Dict1[T[j]]]=1
    end
    return A
end 
A=Adjmat(Dict1,E) #Adjacent Matrix
println(A)

# returns indegree and outdegree of the given vertex
function degree(A,k)
    s = Dict1[k]
    in_ = 0
    out = 0
    for i in 1:length(A[1,:])
        in_ = in_ + A[s,i]
        out = out + A[i,s]
    end
    return in_ , out
end

da =degree(A,'a') #degree of a
db =degree(A,'b') #degree of b
dc =degree(A,'c') #degree of c
dd =degree(A,'d') #degree of d
de =degree(A,'e') #degree of e
println("Degree of a is ", da)
println("Degree of b is ", db)
println("Degree of c is ", dc)
println("Degree of d is ", dd)
println("Degree of e is ", de)

# creates degreee matrix
function degmat(A)
    Deg=[]
    for i in 1:novertex
        X, Y =degree(A,Vertex[i])
        z = X+Y
        append!(Deg,z)
    end
    k = length(Deg)
    B = zeros(Int64,(k,k))
    for j in 1:k
        B[j,j]=Deg[j]
    end
    return B
end

degmat(A) #Degree matrix

Append("ca") #ca
Append("cb") #cb
Append("cd") #cd
Append("ec") #ec
println(E)
A=Adjmat(Dict1,E) #Adjacent Matrix
println(A)
da =degree(A,'a') #degree of a
db =degree(A,'b') #degree of b
dc =degree(A,'c') #degree of c
dd =degree(A,'d') #degree of d
de =degree(A,'e') #degree of e
println("Degree of a is ", da)
println("Degree of b is ", db)
println("Degree of c is ", dc)
println("Degree of d is ", dd)
println("Degree of e is ", de)
degmat(A) #Degree matrix

using GraphRecipes, Plots
graphplot(A, names=1:novertex,arrow=arrow(:head, :closed, 1, 1), curvature_scalar=0.05)


