"""Symbolic source of the three metric forms (brief sections 3 and 12).

X=S^-1, Y=S^-1/2, Z=S^1/2. Physical matrices O are Hermitian; their
coefficient representation XO is S-Hermitian, not generally Euclidean
Hermitian. The anticommutator is an explicitly labelled approximation for
products and spectral functions. SymPy expressions here state the formulas
implemented numerically in models/nonorthogonal.py and bases.py.
"""
from __future__ import annotations

import sympy as sp

H, S, X, Y, Z, V, A, B, q, dH, dS, dH_cell, dS_cell = sp.symbols(
    "H S X Y Z V A B sigma_z dH dS dH_cell dS_cell", commutative=False)
HERMITIAN = {H, S, X, Y, Z, V, A, B, q, dH, dS, dH_cell, dS_cell}


def dagger(expression):
    """Hermitian conjugation: reverse products, conjugate scalars, O^+=O.

    The raw row-centred dipole is excluded: A here is the Hermitian
    midpoint-centred connection returned by siesta.midpoint_connection.
    """
    expression = sp.sympify(expression)
    if expression in HERMITIAN:
        return expression
    if expression.is_commutative:
        return sp.conjugate(expression)
    if expression.is_Add:
        return sp.Add(*(dagger(term) for term in expression.args))
    if expression.is_Mul:
        return sp.Mul(*(dagger(term) for term in reversed(expression.args)))
    if expression.is_Pow and expression.exp.is_Integer:
        return dagger(expression.base)**expression.exp
    return sp.adjoint(expression)


def _reduce_word(word):
    """Reduce adjacent metric factors and commute spin only through the metric."""
    word = list(word)
    rules = {(Y, Y): [X], (Z, Z): [S], (S, X): [], (X, S): [],
             (Y, Z): [], (Z, Y): [], (X, Z): [Y], (Z, X): [Y],
             (S, Y): [Z], (Y, S): [Z]}
    changed = True
    while changed:
        changed = False
        for i in range(len(word)-1):
            pair = tuple(word[i:i+2])
            if pair in rules:
                word[i:i+2] = rules[pair]
                changed = True
                break
            if word[i] == q and word[i+1] in (S, X, Y, Z):
                word[i:i+2] = list(reversed(pair))
                changed = True
                break
    return tuple(word)


def normal_form(expression, *, trace=False):
    """Expand an NC polynomial, reduce YY=X, SX=XS=I and cyclic trace words.

    Z is used solely to express similarity ``f(H_L)=Z f(XH) Y``. The
    rewrite rules also include YZ=ZY=I. No commutation of H, V, A or a
    spectral function with S is assumed. The spin q commutes only with the
    spin-scalar metric, as required by spin_operators.
    """
    expression = sp.expand(expression)
    collected = {}
    for term in sp.Add.make_args(expression):
        coefficient, word = sp.Integer(1), []
        for factor in sp.Mul.make_args(term):
            if factor.is_commutative:
                coefficient *= factor
            elif factor.is_Pow and factor.exp.is_Integer and factor.exp >= 0:
                word.extend([factor.base]*int(factor.exp))
            else:
                word.append(factor)
        reduced = _reduce_word(word)
        if trace and reduced:
            while True:
                rotations = [_reduce_word(reduced[i:]+reduced[:i]) for i in range(len(reduced))]
                best = min(rotations, key=lambda w: (len(w), tuple(map(str, w))))
                if len(best) == len(reduced):
                    reduced = best
                    break
                reduced = best
        collected[reduced] = collected.get(reduced, 0)+coefficient
    return sp.Add(*(sp.simplify(c)*sp.Mul(*word) for word, c in collected.items() if c != 0))


def operator_forms(physical):
    """Return O_c=XO, O_L=YOY, O_sym=(XO+OX)/2 for a physical matrix O."""
    return {"physical": physical, "generalised": X*physical, "lowdin": Y*physical*Y,
            "anticommutator": (X*physical+physical*X)/2}


def hamiltonian_forms():
    """H_c=XH, H_L=YHY, H_sym={X,H}/2; H_sym is not similar to H_c."""
    return operator_forms(H)


def velocity_forms(mode="covariant"):
    """Eq. (4): V=dH-(dS X H+H X dS)/2+i(H X A-A X H).

    bare uses cell derivatives and omits A; berry_connection uses atom
    derivatives and omits A; covariant requires the midpoint dipole A.
    """
    if mode not in ("bare", "berry_connection", "covariant"):
        raise ValueError("unknown velocity mode")
    dh, ds = (dH_cell, dS_cell) if mode == "bare" else (dH, dS)
    physical = dh-(ds*X*H+H*X*ds)/2
    if mode == "covariant":
        physical += sp.I*(H*X*A-A*X*H)
    return operator_forms(physical)


def spin_forms():
    """Sigma_z=S q, Sigma_c=q, Sigma_L=q, Sigma_sym=q when [S,q]=0."""
    return operator_forms(S*q)


def spin_current_forms():
    """J=(V X Sigma+Sigma X V)/2=(V q+q V)/2 for spin-scalar S."""
    return operator_forms((V*q+q*V)/2)


def spin_hall_normalisation():
    """Symbolic source for reconstruct.shc's Pauli-spin Bastin units."""
    bastin = sp.Symbol("KuboBastin", real=True)
    num_cells, cell_area = sp.symbols("N_k A_cell", positive=True)
    return 2*sp.pi*bastin/(num_cells*cell_area)


def hermiticity(forms):
    """Return Euclidean and S-Hermiticity verdicts via symbolic reduction."""
    result = {}
    for name, expression in forms.items():
        result[name] = {
            "hermitian": normal_form(dagger(expression)-expression) == 0,
            "s_hermitian": normal_form(dagger(expression)*S-S*expression) == 0,
        }
    return result


def trace_equivalence():
    """Prove generalised/Lowdin equality for arbitrary catalogue spectral slots.

    H_L=Z (XH) Y and YZ=I imply f(H_L)=Z f(XH) Y by functional calculus
    (polynomial induction, then analytic continuation). Here F_j are
    arbitrary, noncommutative f_j(XH), so this is not a finite-order test.
    Substitute that similarity into the whole word and reduce the cyclic
    polynomial. TD words follow identically; MSD and projectors follow by
    expanding their sums/commutators and applying linearity of the trace.
    """
    result = {"hamiltonian_similarity": normal_form(Z*(X*H)*Y-Y*H*Y) == 0}
    for name, statics, slots in [("dos", 0, 1), ("spectral_operator", 1, 1),
                                  ("greenwood", 2, 2), ("bastin", 2, 2),
                                  ("time_correlation", 2, 3), ("projected_evolution", 3, 3),
                                  ("msd_expanded_term", 2, 3)]:
        operators = sp.symbols(f"O0:{statics}", commutative=False)
        functions = sp.symbols(f"F0:{slots}", commutative=False)
        generalised, lowdin = sp.Integer(1), sp.Integer(1)
        for i in range(max(statics, slots)):
            if i < statics:
                generalised *= X*operators[i]
                lowdin *= Y*operators[i]*Y
            if i < slots:
                generalised *= functions[i]
                lowdin *= Z*functions[i]*Y
        residual = normal_form(lowdin-generalised, trace=True)
        result[name] = {"generalised": generalised, "lowdin_after_similarity": lowdin,
                        "trace_residual": residual, "equal": residual == 0}
    # A single-operator trace is equal. Already with two static operators
    # and no spectral functions the anticommutator is not generally equal.
    single = normal_form(operator_forms(A)["anticommutator"]-X*A, trace=True)
    pair = normal_form(operator_forms(A)["anticommutator"]*operator_forms(B)["anticommutator"]-
                       X*A*X*B, trace=True)
    result["anticommutator"] = {
        "single_operator_equal": single == 0, "pair_residual": pair,
        "generic_products_equal": pair == 0,
        "sufficient_exact_condition": "[S,H]=0 and [S,O_i]=0 for every static slot",
        "classification": "approximation for DOS, Bastin, Greenwood and TD in general",
    }
    return result


def describe():
    """Return printable symbolic forms and verdicts for the one self-check."""
    forms = {"V": velocity_forms(), "Sigma_z": spin_forms(), "J": spin_current_forms()}
    return {name: {"forms": {route: str(normal_form(expr)) for route, expr in item.items()},
                   "hermiticity": hermiticity(item)} for name, item in forms.items()}
