"""
//
// This program contains unit tests for the pauli_manipulations.py file to support code used in:
// Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility.
//
//
"""
# %%
import unittest
import numpy as np
from pauli_manipulations import PauliTerm, PauliH, PauliU

class PauliTermTest(unittest.TestCase):
    """ Unit test cases for PauliTerm """

    ##################################################
    # Init, str, equals, to_pmr_line
    ##################################################
    def test_init_and_str(self):
        x1 = PauliTerm(1.0, [1], ['X'], 1)
        self.assertEqual("X", x1.__str__())
        y2 = PauliTerm(-0.3*1j, [2], ['Y'], 3)
        self.assertEqual("(-0-0.3j)*IYI", y2.__str__())
        xyz = PauliTerm(0.1-0.3*1j, [1,3,4], ['X','Y','Z'], 6)
        self.assertEqual("(0.1-0.3j)*XIYZII", xyz.__str__())

    def test_equals(self):
        # basic 1q tests
        op1 = PauliTerm(1.0, [1], ['X'], 1)
        op2 = PauliTerm(1.0, [1], ['X'], 1)
        self.assertEqual(op1 == op2, True)
        op2 = PauliTerm(0.9, [1], ['X'], 1)
        self.assertEqual(op1 == op2, False)
        op2 = PauliTerm(1.0, [1], ['X'], 2)
        self.assertEqual(op1 == op2, False)
        # 3q tests
        op1 = PauliTerm(1.0, [1, 2], ['X', 'Y'], 3)
        op2 = PauliTerm(1.0, [1, 2], ['X', 'Y'], 3)
        self.assertEqual(op1 == op2, True)
        op2 = PauliTerm(1.0, [1, 3], ['X', 'Y'], 3)
        self.assertEqual(op1 == op2, False)
        op2 = PauliTerm(1.0, [1, 2], ['X', 'Z'], 3)
        self.assertEqual(op1 == op2, False)
        op2 = PauliTerm(1.0, [1, 2], ['X', 'Y'], 2)
        self.assertEqual(op1 == op2, False)

    def test_to_pmr_str(self):
        op = PauliTerm(0.3, [1,3,4], ['X','Y','Z'], 6)
        pmr_str = "0.300000 1 X 3 Y 4 Z"
        self.assertEqual(pmr_str, op.to_pmr_line())

    ##################################################
    # Operations
    ##################################################
    def test_mult(self):
        op = PauliTerm(0.3, [1,3,4], ['X','Y','Z'], 6)
        op = 1.73 * op
        self.assertAlmostEqual(op.c, 1.73 * 0.3)
        op = PauliTerm(-0.6, [1,3,4], ['X','Y','Z'], 6)
        op = 3.14*1j * op
        self.assertAlmostEqual(op.c, 3.14*1j * -0.6)

    def test_dagger(self):
        op = PauliTerm(0.3*1j, [1,3,4], ['X','Y','Z'], 6)
        conj = op.get_dagger()
        self.assertAlmostEqual(conj.c, np.conj(op.c))

    def test_bi_conj(self):
        op1 = PauliTerm(0.3*1j, [1,3,4], ['X','Y','Z'], 6)
        op2 = PauliTerm(1.25, [1, 5], ['Y', 'Z'], 6)
        op3 = PauliTerm(-0.1*1j, [2, 3, 5, 6], ['Y', 'Y', 'X', 'Z'], 6)
        conj_op = op1.bi_conj(op2, op3)
        c = op1.c * op2.c * np.conj(op3.c)
        answer = PauliTerm(c, [1,2,4,5,6], ['Z','Y','Z','Y','Z'], 6)
        self.assertAlmostEqual(conj_op.c, answer.c)
        self.assertEqual(conj_op.supp, answer.supp)
        self.assertEqual(conj_op.ops, answer.ops)
        # next example
        op1 = PauliTerm(0.3*1j, [1,3,4], ['X','Y','Z'], 6)
        op2 = PauliTerm(1.25, [1, 5], ['Y', 'Z'], 6)
        op3 = PauliTerm(-0.1*1j, [2, 3, 6], ['Y', 'Y', 'Z'], 6)
        conj_op = op1.bi_conj(op2, op3)
        c = op1.c * op2.c * np.conj(op3.c)
        answer = PauliTerm(-1j*c, [1,2,4,5,6], ['Z','Y','Z','Z','Z'], 6)
        self.assertAlmostEqual(conj_op.c, answer.c)
        self.assertEqual(conj_op.supp, answer.supp)
        self.assertEqual(conj_op.ops, answer.ops)

class PauliHTest(unittest.TestCase):
    """ Unit test cases for PauliH """

    ##################################################
    # Test init, __str__, and __contains__
    ##################################################
    def test_init_and_str(self):
        # build H = X1 + \lam Z1 model
        lam = 0.25
        x1 = PauliTerm(1.0, [1], ['X'], 3)
        z1 = PauliTerm(1.0, [1], ['Z'], 3)
        h = PauliH(3, [x1, lam * z1])
        h_str = "XII + 0.25*ZII"
        self.assertEqual(h_str, h.__str__())
        # build 2q prl model
        z1z2 = PauliTerm(1.0, [1, 2], ['Z','Z'], 2)
        x1 = PauliTerm(1.0, [1], ['X'], 2)
        x2 = PauliTerm(1.0, [2], ['X'], 2)
        h0 = PauliH(2, [z1z2, 0.1 * x1, 0.1 * x2])
        h0_str = "ZZ + 0.1*XI + 0.1*IX"
        self.assertEqual(h0_str, h0.__str__())
        z1 = PauliTerm(1.0, [1], ['Z'], 2)
        z2 = PauliTerm(1.0, [2], ['Z'], 2)
        h1 = PauliH(2, [z1, z2])
        h1_str = "ZI + IZ"
        self.assertEqual(h1_str, h1.__str__())
        h = h0 + 0.9 * h1
        h_str = h0_str + " + 0.9*ZI + 0.9*IZ"
        self.assertEqual(h_str, h.__str__())

    def test_pmr_str(self):
        # build 2q prl model
        z1z2 = PauliTerm(1.0, [1, 2], ['Z','Z'], 2)
        x1 = PauliTerm(1.0, [1], ['X'], 2)
        x2 = PauliTerm(1.0, [2], ['X'], 2)
        h0 = PauliH(2, [z1z2, 0.1 * x1, 0.1 * x2])
        z1 = PauliTerm(1.0, [1], ['Z'], 2)
        z2 = PauliTerm(1.0, [2], ['Z'], 2)
        h1 = PauliH(2, [z1, z2])
        h = h0 + 0.9 * h1
        pmr_str = "1.000000 1 Z 2 Z\n0.100000 1 X\n0.100000 2 X\n0.900000 1 Z\n0.900000 2 Z"
        self.assertEqual(pmr_str, h.to_pmr_str())

    def test_contains(self):
        z1 = PauliTerm(1.0, [1], ['Z'], 2)
        z2 = PauliTerm(1.0, [2], ['Z'], 2)
        hz = PauliH(2, [z1, z2])
        self.assertEqual(z1 in hz, True)
        self.assertEqual(z2 in hz, True)
        op = PauliTerm(0.9, [1], ['Z'], 2)
        self.assertEqual(op in hz, False)
        x1 = PauliTerm(1.0, [1], ['X'], 2)
        self.assertEqual(x1 in hz, False)
        x2 = PauliTerm(1.0, [2], ['X'], 2)
        hx = PauliH(2, [x1, x2])
        self.assertEqual(x1 in hx, True)
        self.assertEqual(x2 in hx, True)
        xx = PauliTerm(1.0, [1,2], ['X', 'X'], 2)
        self.assertEqual(xx in hx, False)

    ##################################################
    # Test add_term and simplify
    ##################################################
    def test_add_term_and_simplify(self):
        # first suppress simplify with addition
        h = PauliH(3)
        x1 = PauliTerm(1.0, [1], ['X'], 3)
        h.add_term(x1)
        self.assertEqual(x1 in h, True)
        h.add_term(x1, simplify=False)
        self.assertEqual(x1 in h, True)
        h.simplify()
        self.assertEqual(x1 in h, False)
        self.assertEqual(2 * x1 in h, True)
        # now check that simplification occurs with init
        h = PauliH(3, [x1, x1])
        self.assertEqual(x1 in h, False)
        self.assertEqual(2 * x1 in h, True)
        # now check that simplification works with add
        h = PauliH(3)
        h.add_term(x1)
        h.add_term(x1)
        self.assertEqual(x1 in h, False)
        self.assertEqual(2 * x1 in h, True)
        # check simplification with imaginary
        op1 = PauliTerm(0.25+0.1*1j, [1, 3], ['X', 'Y'], 3)
        op2 = PauliTerm(0.3-0.1*1j, [1, 3], ['X', 'Y'], 3)
        op3 = PauliTerm(0.25, [2, 3], ['Z', 'Z'], 3)
        h.add_term(op1, simplify=False)
        h.add_term(op2, False)
        h.add_term(op3, False)
        h.simplify()
        self.assertEqual(op3 in h, True)
        self.assertEqual(op1 in h, False)
        self.assertEqual(op1 in h, False)
        op = PauliTerm(0.25+0.3, [1, 3], ['X', 'Y'], 3)
        self.assertEqual(op in h, True)

class PauliUTest(unittest.TestCase):
    """ Unit test cases for PauliU """

    ##################################################
    # Test init, __str__, and __contains__
    ##################################################
    def test_init_str_contains(self):
        # default should have identity operator
        u = PauliU(3)
        id = PauliTerm(1.0, [], [], 3)
        self.assertEqual(id in u, True)
        self.assertEqual("III", u.__str__())
        # add terms that anti-commute
        c = np.random.normal(0, 1, size=3)
        c /= np.linalg.norm(c)
        op1 = PauliTerm(c[0], [1,2,3], 3*['X'], 3)
        u.add_term(op1)
        op2 = PauliTerm(c[1], [1,2], ['X', 'Z'], 3)
        u.add_term(op2)
        op3 = PauliTerm(c[2], [1], ['Y'], 3)
        u.add_term(op3)
        self.assertEqual(u.validate_normalization(), True)
        self.assertEqual(u.validate_anticommutation(), True)
        self.assertEqual(op1 in u, True)
        self.assertEqual(op2 in u, True)
        self.assertEqual(op3 in u, True)

    ##################################################
    # Test set_as_random and get_inverse
    ##################################################
    def test_set_as_random(self):
        u = PauliU(3)
        u.set_as_random(2*3 + 1)
        ops = ['XXX', 'XXY', 'XXZ', 'XYI', 'XZI', 'YII', 'ZII']
        str_u = u.__str__()
        for op in ops:
            self.assertEqual(op in str_u, True)
        u.set_as_random(3)
        self.assertEqual(u.validate_anticommutation(), True)
        self.assertEqual(u.validate_normalization(), True)
        str_u = u.__str__()
        count = 0
        for op in ops:
            if op in str_u:
                count += 1
        self.assertEqual(count, 3)
        self.assertEqual(u.validate_anticommutation(), True)
        self.assertEqual(u.validate_normalization(), True)

    def test_get_inverse(self):
        u = PauliU(3)
        u.set_as_random(2*3 + 1, eps=np.random.random())
        udg = u.get_inverse()
        for pt in udg.terms:
            conj_pt = PauliTerm(np.conj(pt.c), pt.supp, pt.ops, pt.n)
            self.assertEqual(conj_pt in u, True)


class PauliUConjugateH(unittest.TestCase):
    """ 
    Unit test cases for a PauliU object conjugating
    a PauliH object.
    """

    def test_specific_rotation(self):
        # build 2q prl model
        n = 2
        z1z2 = PauliTerm(1.0, [1, 2], ['Z','Z'], n)
        x1 = PauliTerm(1.0, [1], ['X'], n)
        x2 = PauliTerm(1.0, [2], ['X'], n)
        h0 = PauliH(n, [z1z2, 0.1 * x1, 0.1 * x2])
        z1 = PauliTerm(1.0, [1], ['Z'], n)
        z2 = PauliTerm(1.0, [2], ['Z'], n)
        h1 = PauliH(n, [z1, z2])
        h = h0 + 0.9 * h1
        # rotate it with a hand-crafted unitary
        p1 = PauliTerm(0.3, [1,2], ['X', 'Y'], n)
        p2 = PauliTerm(np.sqrt(1-0.3**2), [1], ['Z'], n)
        u = PauliU(2, [p1, p2])
        uh = h.conjugate(u)
        # form rotated h manually
        a0 = 0.082
        a1 = 0.738
        a2 = 0.515127
        a3 = 0.0572364
        x1 = PauliTerm(-a0, [1], ['X'], n)
        z1 = PauliTerm(a1, [1], ['Z'], n)
        x2 = PauliTerm(a0, [2], ['X'], n)
        yx = PauliTerm(a2, [1,2], ['Y', 'X'], n)
        xy = PauliTerm(a2, [1,2], ['X', 'Y'], n)
        zy = PauliTerm(a3, [1,2], ['Z', 'Y'], n)
        z2 = PauliTerm(a1, [2], ['Z'], n)
        yz = PauliTerm(-a3, [1,2], ['Y', 'Z'], n)
        zz = PauliTerm(1.0, [1,2], ['Z', 'Z'], n)
        man_h = PauliH(n, [x1, z1, x2, yx, xy, zy, z2, yz, zz])
        self.assertEqual(uh == man_h, True)

    def test_self_inverse(self):
        n = 10
        # build 2q prl model
        z1z2 = PauliTerm(1.0, [1, 2], ['Z','Z'], n)
        x1 = PauliTerm(1.0, [1], ['X'], n)
        x2 = PauliTerm(1.0, [2], ['X'], n)
        h0 = PauliH(n, [z1z2, 0.1 * x1, 0.1 * x2])
        z1 = PauliTerm(1.0, [1], ['Z'], n)
        z2 = PauliTerm(1.0, [2], ['Z'], n)
        h1 = PauliH(n, [z1, z2])
        h = h0 + 0.9 * h1
        # generate random u and conjugate h
        u = PauliU(n)
        for _ in range(100):
            l = np.random.choice(range(0, 2*n + 2))
            eps = np.random.random()
            u.set_as_random(l, eps)
            conj_h = h.conjugate(u)
            conj2_h = conj_h.conjugate(u)
            self.assertEqual(conj2_h, h)
# %%

# %%
if __name__ == "__main__":
    unittest.main()
#%%