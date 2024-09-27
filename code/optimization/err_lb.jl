# computes lower bound on the initial guess for the final good firm’s optimization problem

#global variables declaration
global eta::Float64, theta::Float64, deltaa::Float64,  deltay::Float64, psi_ky::Float64, r::Float64, alpha::Float64

function err_lb( ay::Float64, Omega::Float64, rho::Float64, tfp::Float64)
    
    # we declare the variables that are going to be used in the function
    err::Float64
    top::Float64
    bot::Float64        
    cons_p::Float64
    F::Float64
    Fprime::Float64

    # we compute F
    F = exp(-theta*ay^(1.0/theta))

    # we compute F'
    Fprime = -exp(-theta*ay^(1.0/theta)*ay^(1.0/theta-1))

    #we compute err
    err = 1.0 - rho*Omega*(F + Fprime*cons_p)

end
