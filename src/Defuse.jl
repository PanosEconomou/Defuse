module Defuse

greet() = print("Welcome to defuse!")

export folded_minimal_modular_data, build_coefficient_matrix, invert_coefficients, normalize_system

include("folded_modular_data.jl")
include("solve_modular_bootstrap.jl")

end # module Defuse
