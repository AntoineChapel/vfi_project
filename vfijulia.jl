kmin::Int8 = 1
kmax::Int8 = 25
prec::Int16 = 500

kgrid::Vector{Float64} = (collect(kmin:(kmax-kmin)/(prec-1):kmax))
gk::Vector{Float64} = (collect(kmin:(kmax-kmin)/(prec-1):kmax))
Vk0::Vector{Float64} = ones(prec)


const A::Int8 = 10
const alpha::Float64 = 0.5
const beta::Float64 = 0.9

norm::Float64 = 1e5
tol::Float64 = 1e-6
maxiter::Int16 = 1000
n_iter::Int16 = 0 

Vk = Vector{Float64}(undef, prec)
Vkprim = Vector{Float64}(undef, prec)
value_array = Matrix{Float64}(undef, prec, prec)

Vk = Vk0

while n_iter < maxiter && norm > tol
  for iprim in 1:prec
    kprim = kgrid[iprim]
    for i in 1:prec
      k = kgrid[i]
      c = A*(k^alpha) - kprim
      if c > 0 
        value_array[i, iprim] = log(c) + beta*Vk[iprim]
      else
        value_array[i, iprim] = -1e6
      end
    end
  end

  global Vkprim = reduce(max, value_array, dims=2)
  global norm = maximum(abs.(Vkprim - Vk))
  global Vk = Vkprim
  global n_iter += 1 

  println("Iteration: ", n_iter, " norm: ", norm)
  end

  for row in 1:prec
    gk[row] = kgrid[argmax(value_array[row, :])]
  end

kstar = kgrid[argmin(abs.(gk - kgrid))]
println("The steady state value of capital is: ", kstar)
