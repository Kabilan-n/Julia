using LinearAlgebra

function myJacobi(A,x,b,maxIter,errorLimit,resLimit,omega)
    nRow,nCol = size(A)
    row , col = size(transpose(x))
    xold = x
    k = 0
    relError = zeros(maxIter+1,1)
    Notsolved = true
    while Notsolved ==true
        k = k + 1
        xnew = copy(xold)
        
        for i in 1:nRow
            E = xold[i] .+ (b[i] .- reshape(A[i,:],row, col)*xold)/ A[i,i]
            xnew[i] = E[1]
        end
        
        currentError = norm(xnew - xold)
        relError[k] = currentError / norm(xnew)
        
        if norm(b - A*xnew) <= resLimit || currentError <= errorLimit  || k > maxIter
            Notsolved = false
        else
            xold = xnew
        end
    end
    
    return xnew,k,relError
end

A = [-4 2 1 0 0;1 -4 1 1 0;2 1 -4 1 2; 0 1 1 -4 1;0 0 1 2 -4]
A = float(A)

b = [-4 11 -16 11 -4]
b = float(transpose(b))

maxIter = 30
errorLimit = 0.00001
resLimit = 0.00001
omega =1

x = [1.0 1.0 1.0 1.0 1.0]
x = float(x)
x =transpose(x)

xnew,k,relError = myJacobi(A,x,b,maxIter,errorLimit,resLimit,omega)
println("Number of iterations")
println(k)
println("the solution vector is")
println(transpose(xnew))
println("recomputed b is")
println(transpose(A*xnew))
println("original b is")
println(b)

function mySOR(A,x,b,maxIter,errorLimit,resLimit,omega)
    nRow,nCol = size(A)
    xold = x
    row , col = size(transpose(x))
    k = 0
    relError = zeros(maxIter+1,1)
    Notsolved = true
    
    while Notsolved ==true
        k = k + 1
        xnew = xold
        
        for i in 1:nRow
            E = xnew[i] .+ omega*(b[i] .- reshape(A[i,:],row,col)*xnew)/ A[i,i]
            xnew[i] = E[1]
        end
        
        currentError = norm(xnew - xold)
        relError[k] = currentError / norm(xnew)
        
        if norm(b - A*xnew) <= resLimit || currentError <= errorLimit  || k > maxIter
            Notsolved = false
        else
            xold = xnew
        end
    end
    
    return xnew,k,relError
end

A = [-4 2 1 0 0;1 -4 1 1 0;2 1 -4 1 2; 0 1 1 -4 1;0 0 1 2 -4]
A = float(A)

b = [-4 11 -16 11 -4]
b = float(transpose(b))

maxIter = 30
errorLimit = 0.00001
resLimit = 0.00001
omega =1.6

x = [1.0 1.0 1.0 1.0 1.0]
x = float(x)
x =transpose(x)

xnew,k,relError = mySOR(A,x,b,maxIter,errorLimit,resLimit,omega)
println("Number of iterations")
println(k)
println("the solution vector is")
println(transpose(xnew))
println("recomputed b is")
println(transpose(A*xnew))
println("original b is")
println(b)


