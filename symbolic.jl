
using RigidBodyDynamics
using StaticArrays
using LinearAlgebra: cross 
using SymPy

## Walker parameters
a = 0.5
b = 0.5
m₁ = 0.5 
m₂ = 0.5
mₕ = 10
g = 9.8
@syms α t q̇₁₋ q̇₁₊ q̇₂₋ q̇₂₊ real = true
@symfuns q₁ q₂ q̇₁ q̇₂ q̈₁ q̈₂
## Swing leg phase, figure 1
TF_1_wrt_0 = @SMatrix[cos(-q₁(t))      sin(-q₁(t))       b*sin(-q₁(t))       ;
                      -sin(-q₁(t))     cos(-q₁(t))       b*cos(-q₁(t))       ;
                      0                0                 1                   ]
TF_h_wrt_1 = @SMatrix[1                0                 0                   ;
                      0                1                 a                   ;
                      0                0                 1                   ]
TF_2_wrt_h = @SMatrix[cos(q₂(t)-q₁(t)) -sin(q₂(t)-q₁(t)) a*sin(q₂(t)-q₁(t))  ;
                      sin(q₂(t)-q₁(t)) cos(q₂(t)-q₁(t))  -a*cos(q₂(t)-q₁(t)) ;
                      0                0                 1                   ]
TF_1_wrt_0 = simplify.(TF_1_wrt_0)
TF_h_wrt_0 = simplify.(TF_1_wrt_0 * TF_h_wrt_1)
TF_2_wrt_0 = simplify.(TF_h_wrt_0 * TF_2_wrt_h)
p₁ = TF_1_wrt_0[7:8]
pₕ = TF_h_wrt_0[7:8]
p₂ = TF_2_wrt_0[7:8]
v₁ = diff.(p₁,t)
vₕ = diff.(pₕ,t)
v₂ = diff.(p₂,t)
K = simplify(0.5*(m₁*sum(v₁.^2) + m₂*sum(v₂.^2) + mₕ*sum(vₕ.^2)))
P = simplify((m₁*p₁[2]+ m₂*p₂[2] + mₕ*pₕ[2])*g)
## Lagrangian = simplify(K-P)
function manipulator_matrices(K::Sym, P::Sym, q::Array{Sym,1})
    n = length(q)
    m = Array{Sym}(undef, n, n)
    c = Array{Sym}(undef, n, n)
    g = Array{Sym}(undef, n)
    for i = 1:n
        for j in 1:n
            m[i,j] = diff(diff(K,sympy.Derivative(q[i], t)),
                          sympy.Derivative(q[j], t))
        end
    end
    m = simplify.(m)
    for i = 1:n
        for j = 1:n
            sum = 0
            for k = 1:n
                sum += diff(m[i,j],q[k])
                sum += diff(m[i,k],q[j])
                sum -= diff(m[k,j],q[i])
                sum *= sympy.Derivative(q[k], t)
            end
            c[i,j] = 0.5*sum
        end
    end
    c = simplify.(c)
    for i in 1:n
        g[i] = diff(P,q[i])
    end
    g = simplify.(g)
    return m,c,g
end
M,C,G = manipulator_matrices(K,P,[q₁(t),q₂(t)])
## Impact phase, figure 2
p₂₊ = [-b*sin(q₂(t));
       b*cos(q₂(t)) ;
       0            ]
pₕ₊ = [-(a+b)*sin(q₂(t));
       (a+b)*cos(q₂(t)) ;
       0                ]
p₁₊ = pₕ₊ +[-a*sin(-q₁(t)) ;
            -a*cos(-q₁(t)) ;
            0              ]
if(length(p₁) < 3)
    push!(p₁,0)
    push!(pₕ ,0)
    push!(p₂,0)
end
p₁₋ = p₁ 
pₕ₋ = pₕ 
p₂₋ = p₂
v₁₋ = cross([0; 0; q̇₁₋], p₁₋)
vₕ₋ = cross([0; 0; q̇₁₋], pₕ₋)
v₂₋ = vₕ₋ + cross([0; 0; q̇₂₋], p₂₋ - pₕ₋)
vₕ₊ = cross([0; 0; q̇₂₊], pₕ₊)
v₂₊ = cross([0; 0; q̇₂₊], p₂₊)
v₁₊ = vₕ₊ + cross([0; 0; q̇₁₊], p₁₊ - pₕ₊)

## conservation of angular momentum
conserved₀ = simplify.(cross(p₁₋,m₁*v₁₋) + cross(p₂₋,m₂*v₂₋) +
                       cross(pₕ₋,mₕ*vₕ₋) + ( cross(p₁₊,m₁*v₁₊) +
                                             cross(p₂₊,m₂*v₂₊) +
                                             cross(pₕ₊,mₕ*vₕ₊)))[3]
p₁ₕ = p₁ - pₕ
conservedₕ = simplify.(cross(p₁ₕ,m₁*v₁₋) + cross(p₁ₕ,m₁*v₁₊))[3]
Q₊ = [conserved₀.coeff(q̇₁₊) conserved₀.coeff(q̇₂₊);
      conservedₕ.coeff(q̇₁₊) conservedₕ.coeff(q̇₂₊)]
Q₋ = [conserved₀.coeff(q̇₁₋) conserved₀.coeff(q̇₂₋);
      conservedₕ.coeff(q̇₁₋) conservedₕ.coeff(q̇₂₋)]
simplify.(inv(Q₊)*Q₋)
