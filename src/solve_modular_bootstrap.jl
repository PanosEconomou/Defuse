# takes s and z matrices and spits out the traces of simple lines on primaries
using LinearAlgebra

function build_coefficient_matrix(z,s)
    st              = s'
    nonzero         = findall(!iszero,z)
    coefficients    = Matrix{ComplexF64}(undef, size(z,1)^2, length(nonzero))

    for (i, pos) in enumerate(nonzero)
        coefficients[:,i]  = vec(st[:,pos[1]] * transpose(s[pos[2],:]))
    end

    return coeff
end

function invert_coefficients(coefficients)
    # Get the pivot rows
    qr_decomposition    = qr(transpose(coefficients), ColumnNorm())
    coefficent_rank     = rank(coefficients, atol=1e-10)
    pivot_rows          = qr_decomposition.p[1:coefficent_rank]

    # Invert the matrix of pivot rows
    return inv(coefficients[pivot_rows,:])
end

function normalize_system(coefficients, digits::Int = 14)
    inverse     = invert_coefficients(coefficients)
    system      = real(round.(coefficients * inverse, digits=digits))
    
    # Remove duplicates
    system      = system[unique(i -> system[i,:], 1:size(system,1)), :]
end

