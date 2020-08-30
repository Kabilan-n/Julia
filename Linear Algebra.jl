x = [ -1.1, 0.0, 3.6, -7.2 ]

length(x)

y = [ -1.1; 0.0; 3.6; -7.2 ]

length(y)

a = [1,2]

b = (1,2)

x = [ -1.1, 0.0, 3.6, -7.2 ];

x = [ -1.1, 0.0, 3.6, -7.2 ]

x[3]

x[3] = 14

x

y = x

z = copy(x)

x = [1]

x = 1

x1 = [ 1, -2 ];  x2 = [ 1, 1 ];
X = [x1; x2]

X = [x1, x2]

length(X)

Y = vcat(x1, x2)

Z = hcat(x1, x2)

x = [ 9, 4, 3, 0, 5 ]

x[2:5]

x[4:5] = [ -2, -3 ];

x

x = [ 1, 0, 0, -2, 2 ];
d = x[2:end] - x[1:end-1]

x[2:end]

x[1:end-1]

x = [1, 0];
y = [1, -1];
z = [0, 1];

lst = [x, y, z]

lst[2]

zeros(3)

zeros(3,3)

ones(4)

ones(3,3)

rand(2)

rand(2,2)

randn(2)

randn(2,2)

rand(Int,2)

rand(Int,2,3)

rand(1:10)

rand(1:10,3,3)

import Pkg; 
Pkg.add("Plots")

using Plots

x = rand(1:10,10,1);
plot(x, marker = :circle, legend = false, grid = true)

savefig("x.pdf")

[ 0, 7, 3 ] + [ 1, 2, 0 ]

[ 1, 9 ] - [ 1, 1]

x = [ 0, 2, -1 ]
2 * x

x * 2

x / 3

3 \ x

x - 3

x .- 3

1.5 + [1 -1]

1.5 .+ [1 -1]

x1 = [2 4 6];
x2 = [1 2 3];

x1 / x2

x1 ./ x2

x1 = [1 2 3];
x2 = [1 2 2];

x1 == x2

x1 .== x2

x = [1 2 -1];
abs.(x)

a = [ 1, 2 ]; b = [ 3, 4 ];

alpha = -0.5; beta = 1.5;

c = alpha*a + beta*b

function lincomb(coeff, vectors)
    n = length(vectors[1])
    a = zeros(n);
    for i = 1:length(vectors)
        a = a + coeff[i] * vectors[i];
    end
    return a
end

lincomb( ( -0.5, 1.5), ( [1, 2], [ 3, 4]) )

function lincomb(coeff, vectors)
    return sum( coeff[i] * vectors[i] for i = 1:length(vectors) )
end

lincomb( ( -0.5, 1.5), ( [1, 2], [ 3, 4]) )

x = [ -1, 2, 2 ];
y = [ 1, 0, -3 ];

x'*y

a = rand(2,2);
a^2

a = randn(10^5);
@time a'*a

f(x) = x[1] + x[2] - + x[3] - x[4]^2
f([-1,0,1,2])

M = [1 2 2; 1 2 2; 2 3 4]

nullspace(M)

A = [1 3  1; 1 1 -1; 3 11 6];
b = [9, 1, 35];

m = [148.73, -18.85]; 
c = 54.40;

y_hat(x) = x'*m + c;

x = [0.846, 1]; 
y = 115;

y_hat(x), y

x = [ 2, -1, 2 ];

using LinearAlgebra
norm(x)

norm(x,1)

norm(x,2)

norm(x,Inf)

cos(30)

sin(30)

A = [1 2 3 4;4 5 6 7;7 8 9 10]

m, n = size(A)

reshape(A, (4,3))

H = [0 1 -2; 1 2 -1; 3 0 6]

diag(H)

A = [0 2 -1; -2 1 1]
x = [2, 1, -1]

A*x

a = [1,1]; 
b = [2,-1,1];

import Pkg; 
Pkg.add("DSP")

using DSP

conv(a,b)

a = [1,1];
b = [2,-1,1]; 
c = [1,1,-2];

d = conv(conv(a,b),c)

A = [-1.5 3 2; 1 -1 0];
B = [-1 -1; 0 -2; 1 0];

C = A*B

Q, R = qr(A);

Q

R

A = [1 -2 3; 0 2 2; -4 -4 -4]

B = inv(A)

C = pinv(A)

R = triu(randn(4,4))

R = tril(randn(4,4))

A = randn(10,10); 
b = randn(10);

x1 = A\b

x2 = inv(A)*b

A = rand(1:10,3,3)

rank(A)

svd(A)

lu(A)

eigvals(A)

eigvecs(A)

using PyPlot
x = rand(8,8)
pcolor(x)

using Pkg
Pkg.add("Convex")
Pkg.add("SCS")

using Convex, SCS

# Generate random problem data
m = 4;  n = 5
A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n)

# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# This can be done by: minimize(objective, constraints)
problem = minimize(sumsquares(A * x - b), [x >= 0])

# Solve the problem by calling solve!
solve!(problem, SCS.Optimizer)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimal value
problem.optval

Pkg.add("Distributions")

# Generate data.
n = 2; # dimensionality of data
C = 10; # inverse regularization parameter in the objective
N = 10; # number of positive examples
M = 10; # number of negative examples

using Distributions: MvNormal
# positive data points
pos_data = rand(MvNormal([1.0, 2.0], 1.0), N);
# negative data points
neg_data = rand(MvNormal([-1.0, 2.0], 1.0), M);

function svm(pos_data, neg_data, solver=() -> SCS.Optimizer(verbose=0))
    # Create variables for the separating hyperplane w'*x = b.
    w = Variable(n)
    b = Variable()
    # Form the objective.
    obj = sumsquares(w) + C*sum(max(1+b-w'*pos_data, 0)) + C*sum(max(1-b+w'*neg_data, 0))
    # Form and solve problem.
    problem = minimize(obj)
    solve!(problem, solver)
    return evaluate(w), evaluate(b)
end;

w, b = svm(pos_data, neg_data);

w

b


