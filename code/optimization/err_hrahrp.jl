# Purpose: Solve the rental-houseing firm's optimization problem from the first orderconditions

# # global variables declariation
global eta::Float64, theta::Float64, deltah::Float64, deltaa::Float64, psi_kr::Float64, r::Float64, alpha::Float64, Ah::Float64

function err_hrahrp(ar::Float64, omega::Float64, rho::Float64, tfp::Float64)
    # need to pass tfp so that you can call same rtbis function

    # we declare variables of type doubles
    top::Float64 = 0.0
    bot::Float64 = 0.0
    ph::Float64 = 0.0
    F::Float64 = 0.0
    Fprime::Float64 = 0.0
    cons_p::Float64 = 0.0
    cons_a::Float64 = 0.0
    err::Float64 = 0.0
    fochra::Float64 = 0.0

    # Calculate F using the given formula
    F = exp(-theta * ar^(1.0 / theta))
    # Calculate the derivative of F with respect to ar
    Fprime = -exp(-theta * ar^(1.0 / theta)) * ar^(1.0 / theta - 1.0)
    cons_p = -ar # Consumption of p is negative ar
    cons_a = 1.0 # Consumption of a is constant 1.0

    # ph
    top = r+deltah + rho*omega*(1.0-deltah - psi_kr)*(F+Fprime*cons_p)
    bot = Ah*(1.0 - rho*omega*(F + Fprime*cons_p))
    ph = top/bot

    # foc for adaptative rental housing
    fochra = -ro*(1.0+eta)*(ph*Ah + 1.O - delta -psi_kr)*omega*Fprime*cons_a - r - deltaa

    err = fochra

    return err
end



