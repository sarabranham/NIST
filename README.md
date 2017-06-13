# NIST
# EQNS:
#- transverse components of e/m fields
#    E(t) = [ Ce^(j(ωt - βz)) + De^(j(ωt + βz)) ] * F(t)                  >>>   [C + D] * F(t)
#    H(t) = [ Ce^(j(ωt - βz)) - De^(j(ωt + βz)) ] * Y * k x F(t)          >>>   [C - D] * Y * k x F(t)
#
#    where
#       c/d = complex constants                         F(t) = perp to k, how E(t) varies over cross section
#       k = unit vector in z (axis of transmission)     k x F = rotate F in transverse plane (so H perp E)
#       Y = function of dielectric/mag props.           time dependence = e^jωt and z dep. = e^±jβz
#
#- velocity of wave
#    v = Δz/Δt = ω/β
#    β = 2πf/v = 2π/λ
#
#- power normalization (derivation makes no sense)
#    P = K(|A|^2 - |B|^2)
#
#    if low frequency
#       P = (|A|^2 - |B|^2) / Z.
#
#- voltage and current
#    V = A + B and I = (A-B)/Z.                          where z. = z @ t = 0
#    A = 1/2(V + IZ.) and B = 1/2(V - IZ.)
#
#- boundary conditions
#    E = [ Ae^(j(ωt - βz)) + Be^(j(ωt + βz)) ]
#
#    if a = b (standing wave, min bz + θ = nπ, max bz + θ = (2n+1)π/2
#        E = [ ae^(-jβz) + ce^(jβz) ]           >>>>> MATH IT >>>>>          E = 2ae^(jθ)cos(βz + θ)
#
#- reflection
#    a = bΓ
#    where
#        a = incidence amplitude, b = emerging amplitude, gamma = reflection coefficient
#
#    p(abs) = |b^2| - |a^2|     OR p(abs) = |b^2|(1 - |Γ|^2)
#
#    gamma = 0 either means totally absorbed by termination or something else went wrong
#    gamma = 1 means totally reflected, negative numbers depending on where you define 0