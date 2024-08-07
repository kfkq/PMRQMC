"""
//
// This program contains code generate and apply the "special random unitaries" described in Appendix A of:
// Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility.
//
//
"""
# %%
import numpy as np

class PauliTerm:
    def __init__(self, c, supp, ops, n):
        if len(supp) != len(ops):
            raise ValueError("The length of 'supp' must match the length of 'ops'")
        self.c = c
        # sort the support from smallest to larget always
        supp = np.array(supp)
        ops = np.array(ops)
        sort_idx = np.argsort(supp)
        self.supp = tuple(supp[sort_idx])
        self.ops = tuple(ops[sort_idx])
        self.n = n

    def __str__(self):
        full_ops = ['I'] * self.n
        for idx, op in zip(self.supp, self.ops):
            full_ops[idx - 1] = op
        ops_str = ''.join(full_ops)
        return f"{self.c}*{ops_str}" if self.c != 1 else f"{ops_str}"

    def __eq__(self, other):
        n_check = (self.n == other.n)
        c_check = np.isclose(self.c, other.c)
        sup_check = (self.supp == other.supp)
        op_check = (self.ops == other.ops)
        return n_check and c_check and sup_check and op_check

    def to_pmr_line(self):
        parts = [f"{self.c:.6f}"]
        for idx, op in zip(self.supp, self.ops):
            parts.append(str(idx))
            parts.append(op)
        return ' '.join(parts)

    def __mul__(self, m):
        """Handle scalar multiplication from the left: other * self."""
        if isinstance(m, (int, float, complex)):  # Ensure other is a number
            return PauliTerm(m * self.c, self.supp, self.ops, self.n)
        return NotImplemented

    def __rmul__(self, m):
        """Handle scalar multiplication from the right: self * other."""
        return self.__mul__(m)  # Multiplication is commutative for scalars

    def get_dagger(self):
        """
        Returns p^{\dagger}.
        """
        return PauliTerm(np.conj(self.c), self.supp, self.ops, self.n)

    def bi_conj(self, p2, p3):
        """Conjugate self (p1) by p2 from left and p3^dagger from right to form p4 = p2 * p1 * p3^dagger"""
        if self.n != p2.n or self.n != p3.n:
            raise ValueError("All Pauli terms must have the same number of qubits.")

        # Initialize new coefficient
        new_c = self.c * p2.c * np.conj(p3.c)  # Include conjugate for p3

        # Initialize the result ops
        new_ops = ['I'] * self.n

        # Apply the bi-conjugation across all qubits
        for i in range(self.n):
            # exact operators
            if i + 1 in self.supp:
                op1 = self.ops[self.supp.index(i + 1)]
            else:
                op1 = 'I'
            if i + 1 in p2.supp:
                op2 = p2.ops[p2.supp.index(i + 1)]
            else:
                op2 = 'I'
            if i + 1 in p3.supp:
                op3 = p3.ops[p3.supp.index(i + 1)]
            else:
                op3 = 'I'
            
            # Apply p2 * p1
            interim_op, interim_phase = self.pauli_single_multiply(op2, op1)
            # Apply result * p3^dagger
            new_op, phase = self.pauli_single_multiply(interim_op, op3)
            new_ops[i] = new_op
            new_c *= interim_phase * phase  # Adjust coefficient with combined phases

        new_supp = [i + 1 for i, op in enumerate(new_ops) if op != 'I']
        non_id_ops = [p for p in new_ops if p != 'I']

        return PauliTerm(new_c, new_supp, non_id_ops, self.n)


    @staticmethod
    def pauli_single_multiply(op1, op2):
        """Multiply two single-qubit Pauli operators and return the result and phase."""
        if op1 == 'I':
            return op2, 1
        if op2 == 'I':
            return op1, 1
        if op1 == op2:
            return 'I', 1  # Multiplying same operators results in the identity

        # Mapping rules for multiplying Pauli matrices, considering phases
        rules = {
            ('X', 'Y'): ('Z', 1j),
            ('Y', 'Z'): ('X', 1j),
            ('Z', 'X'): ('Y', 1j),
            ('Y', 'X'): ('Z', -1j),
            ('Z', 'Y'): ('X', -1j),
            ('X', 'Z'): ('Y', -1j),
        }
        result, phase = rules.get((op1, op2), rules.get((op2, op1), ('I', -1)))  # Check both tuple orders
        return result, phase

class PauliH:
    def __init__(self, n=None, terms=None):
        self.terms = []
        self.n = n
        if terms is not None:
            for term in terms:
                self.add_term(term, simplify=False)
            self.simplify()

    def __contains__(self, pauli):
        for op in self.terms:
            if op == pauli:
                return True
        return False

    def __eq__(self, other):
        # simplify first to ensure things are correct
        self.simplify()
        other.simplify()
        if len(other.terms) != len(self.terms):
            return False
        for o in other.terms:
            if o not in self:
                return False
        return True

    def add_term(self, term, simplify=True):
        if self.n is not None and term.n != self.n:
            raise ValueError("The number of qubits in the term must match the Hamiltonian.")
        if self.n is None:
            self.n = term.n
        self.terms.append(term)
        if simplify:
            self.simplify()

    def simplify(self):
        """Combine like terms and ensure the coefficients are real at the end."""
        combined_terms = {}
        for term in self.terms:
            key = (term.supp, term.ops)
            if key in combined_terms:
                combined_terms[key].c += term.c
                combined_terms[key].c = self.truncate_imaginary(combined_terms[key].c)
            else:
                combined_terms[key] = PauliTerm(term.c, term.supp, term.ops, self.n)

        new_terms = []
        for term in combined_terms.values():
            if not np.isclose(term.c.imag, 0):
                raise ValueError("The final coefficients must be real.")
            elif not np.isclose(term.c.real, 0):
                new_terms.append(PauliTerm(term.c.real, term.supp, term.ops, self.n))

        self.terms = new_terms

    def __str__(self):
        """Return a string representation of the Hamiltonian."""
        if not self.terms:
            return "0"
        return " + ".join(str(term) for term in self.terms)

    def to_pmr_str(self):
        """Return the PMR formatted string of the Hamiltonian, each term on a new line."""
        return "\n".join(term.to_pmr_line() for term in self.terms)

    def set_as_1q_model(self, diag_pert, lam):
        self.terms = []
        if diag_pert is True:
            self.add_term(PauliTerm(1.0, [1], ['X'], self.n))
            self.add_term(PauliTerm(lam, [1], ['Z'], self.n))
        else:
            self.add_term(PauliTerm(lam, [1], ['X'], self.n))
            self.add_term(PauliTerm(1.0, [1], ['Z'], self.n))

        return 

    def conjugate(self, u):
        if not isinstance(u, PauliU):
            raise TypeError("The object must be an instance of PauliU.")
        new_hamiltonian = PauliH(n=self.n)
        for term_h in self.terms:
            for term_u1 in u.terms:
                for term_u2 in u.terms:
                    # Conjugate h by u1 and u2^dagger
                    new_term = term_h.bi_conj(term_u1, term_u2)
                    new_hamiltonian.add_term(new_term, simplify=False)
        new_hamiltonian.simplify()
        return new_hamiltonian

    def __add__(self, other):
        if not isinstance(other, PauliH):
            return NotImplemented
        if self.n != other.n:
            raise ValueError("Cannot add Hamiltonians with different numbers of qubits.")
        new_terms = self.terms + other.terms
        return PauliH(self.n, new_terms)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            new_terms = [other * term for term in self.terms]
            return PauliH(self.n, new_terms)
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    @staticmethod
    def truncate_imaginary(c, atol=1e-9):
        if np.isclose(c.imag, 0, atol=atol):
            return c.real
        return c

class PauliU:
    def __init__(self,  n=None, terms=None):
        self.terms = []
        self.n = n
        if terms is not None:
            for term in terms:
                self.add_term(term)
            self.validate_normalization()
        if terms is None and n is not None:
            self.terms = [PauliTerm(1.0, [], [], n)]

    def __contains__(self, pauli):
        for op in self.terms:
            if op == pauli:
                return True
        return False

    def add_term(self, term):
        if self.n is not None and term.n != self.n:
            raise ValueError("The number of qubits in the term must match the Hamiltonian.")
        if self.n is None:
            self.n = term.n
        # remove identity term if present
        id = PauliTerm(1.0, [], [], self.n)
        if id in self:
            self.terms = []
        self.terms.append(term)
        self.validate_anticommutation()

    def __str__(self):
        if not self.terms:
            return "0"
        return " + ".join(str(term) for term in self.terms)

    def to_pmr_str(self):
        return "\n".join(term.to_pmr_line() for term in self.terms)

    def validate_normalization(self):
        """Ensure the sum of squares of the coefficients is 1."""
        c = np.array([p.c for p in self.terms])
        if not np.isclose(np.linalg.norm(c), 1):
            raise ValueError("The coefficients are not normalized properly.")
        return True

    def validate_anticommutation(self):
        """Ensure all Pauli terms anti-commute."""
        for i in range(len(self.terms)):
            for j in range(i + 1, len(self.terms)):
                if not self.check_anticommutation(self.terms[i], self.terms[j]):
                    raise ValueError("Not all Pauli terms anti-commute.")
        return True

    def set_as_random(self, l, eps=0.05, mu=0, sig=1, seed=None):
        """
        Forms random n qubit unitary consisting of
        sum of l Pauli, anticommuting Paulis.
        See [J. Chem. Theory Comput. 2020, 16, 1, 190–195].
        """
        if seed is not None:
            np.random.seed(seed)
        if l > 2 * self.n + 1:
            raise ValueError("l cannot be bigger than 2n+1")
        elif l == 0:
            print("No setting done since l = 0")
            return
        self.terms = []
        # choose l random anti-commuting paulis
        canon_paulis = self.form_canon_anticomm_paulis(self.n)
        if l == 1:
            l_paulis = [canon_paulis[0]]
        elif l == 2:
            l_paulis = [canon_paulis[0], canon_paulis[-1]]
        else:
            l_pauli_idx = np.random.choice(range(1, len(canon_paulis)-1), l-2, replace=False)
            l_paulis = [canon_paulis[0], canon_paulis[-1]] + [canon_paulis[i] for i in l_pauli_idx]
        # choose l real, normalized coefficients
        if l == 1:
            c = [1.0]
        else:
            if eps is not None:
                c1 = [np.sqrt(1 - eps)]
                c_rest = np.random.normal(mu, sig, l-1)
                fac = np.sqrt(eps) / np.linalg.norm(c_rest)
                c = c1 + list(fac * c_rest)
            else:
                c = np.random.normal(mu, sig, l)
        c = np.array(c)
        c /= np.linalg.norm(c)
        # build and add as PauliTerms
        for k in range(l):
            p = self.form_pauli_term(c[k], l_paulis[k])
            self.add_term(p)
        self.validate_normalization()
        self.validate_anticommutation()

    def get_inverse(self):
        """
        Returns u^{\dagger}, i.e., inverse of [self].
        """
        new_u = PauliU(self.n, [p.get_dagger() for p in self.terms])
        return new_u

    @staticmethod
    def check_anticommutation(term1, term2):
        differing_positions = 0
        for i in range(term1.n):
            op1 = term1.ops[term1.supp.index(i+1)] if i+1 in term1.supp else 'I'
            op2 = term2.ops[term2.supp.index(i+1)] if i+1 in term2.supp else 'I'
            if op1 != 'I' and op2 != 'I' and op1 != op2:
                differing_positions += 1
        return differing_positions % 2 == 1

    @staticmethod
    def form_canon_anticomm_paulis(n):
        """
        See 4.2 in [https://arxiv.org/pdf/1909.08123.pdf],
        construction (1).
        """
        g = ["X", "Y", "Z"]
        for _ in range(n-1):
            x_part = ["X" + p for p in g]
            y_part = ["Y" + "I" * len(g[0])]
            z_part = ["Z" + "I" * len(g[0])]
            g = x_part + y_part + z_part

        return g

    @staticmethod
    def form_pauli_term(c, pstr):
        """
        Forms PauliTerm objct given
        coefficient [c] and Pauli string
        [pstr].
        """
        supp = []
        ops = []
        for k in range(len(pstr)):
            if pstr[k] != "I":
                supp.append(k + 1)
                ops.append(pstr[k])

        return PauliTerm(c, supp, ops, len(pstr))
# %%