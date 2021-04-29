import sympy as sym
import numpy as np
import itertools


class FuzzyIntegrals:

    def __init__(self, mu=None, f=None):
        self.mu = None
        self.f = None
        self.f_sorted = None
        self.f_sorted_rev = None
        self.lmbd = None
        self.g_X = {}
        self.g_A = []

        if mu is not None and f is not None:
            self.set_params(mu, f)

    def set_params(self, mu, f):
        self.mu = {k: v for k, v in enumerate(mu)}
        self.f = {k: v for k, v in enumerate(f)}
        self.f_sorted = {k: v for k, v in sorted(self.f.items(), key=lambda item: item[1], reverse=False)}
        self.f_sorted_rev = {k: v for k, v in sorted(self.f.items(), key=lambda item: item[1], reverse=True)}
        self.calc_lambda()
        self.calc_g_A()

    def validate_params(self, fn, mu, f):
        if mu is not None and f is not None:
            self.set_params(mu, f)
        else:
            if self.mu is None or self.f is None:
                print("[!] Error: Necessary Inputs Not Provided (mu, f)\t\t--Aborting")
                exit(-1)
            else:
                if mu is not None or f is not None:
                    print(f"[!] Warning: Incomplete Inputs Provided in Calling {fn} Function (mu, f)\t\t--Ignoring")

    def calc_lambda(self):
        eq = ""
        lm = []

        # BUILD LAMBDA EQUATION-----------
        for i in self.mu:
            eq += f"(1+{self.mu[i]}*x)"
            if i != len(self.mu) - 1:
                eq += "*"
        eq += "-x-1"

        # PARSE AND SOLVE LAMBDA EQUATION-----------
        parsed = sym.parse_expr(eq)
        x = sym.Symbol('x')
        lm_eq = sym.Eq(parsed, 0)
        lm_roots = sym.solve(lm_eq, [x])

        # FIND REAL ROOTS ------------------
        for r in lm_roots:
            if sym.sympify(r).is_real:
                lm.append(r)

        # CHECK VALIDITY OF LAMBDA VALUE ----------------------
        lmbd = None
        if np.sum(list(self.mu.values())) > 1:
            for r in lm:
                if -1 < r < 0:
                    lmbd = r
                    break
        elif np.sum(list(self.mu.values())) < 1:
            for r in lm:
                if r > 0:
                    lmbd = r
                    break
        else:
            lmbd = 0
        self.lmbd = lmbd
        return lmbd

    def calc_g_A(self):

        # CALCULATE g(A_i) -------------------------
        states = list(range(len(self.mu)))
        subsets = []
        for s in range(1, len(states) + 1):
            subsets += itertools.combinations(states, s)

        self.g_X = {}
        for s in subsets:
            if len(s) == 1:
                self.g_X[s] = self.mu[s[0]]
            else:
                t = self.combinator(self.mu[s[0]], self.mu[s[1]], self.lmbd)
                for u in s[2:]:
                    t = self.combinator(t, self.mu[u], self.lmbd)
                self.g_X[s] = t

        order = []
        for i, k in enumerate(self.f_sorted_rev):
            fs_copy = self.f_sorted_rev.copy()
            tp = (k,)
            fs_copy.pop(k)
            for j, kk in enumerate(fs_copy):
                if self.f_sorted_rev[k] == self.f_sorted_rev[kk]:
                    tp += (kk,)
            order.append(tuple(sorted(tp)))

        self.g_A = []
        for i, o in enumerate(order):
            tt = o
            for j in range(i):
                tt += order[j]
            tt = tuple(set(tt))
            self.g_A.append(self.g_X[tt])
        self.g_A.reverse()

    def sugeno(self, mu=None, f=None):
        self.validate_params("Sugeno", mu, f)

        # CALCULATE SUGENO INTEGRAL ------------------
        lst = []
        for i, k in enumerate(self.f_sorted):
            lst.append(min(self.f_sorted[k], self.g_A[i]))
        S = max(lst)
        return S

    def choquet(self, mu=None, f=None):
        self.validate_params("Choquet", mu, f)

        # CALCULATE CHOQUET INTEGRAL ------------------
        lst = []
        f_keys = list(self.f_sorted.keys())
        for i, k in enumerate(self.f_sorted):
            if i == 0:
                lst.append((self.f_sorted[f_keys[i]]-0) * self.g_A[i])
            else:
                lst.append((self.f_sorted[f_keys[i]] - self.f_sorted[f_keys[i-1]]) * self.g_A[i])
        C = sum(lst)
        return C

    def get_params(self):
        if self.lmbd is None or not bool(self.g_X):
            print("[!] Error: Object Not Initialized (mu, f)\t\t--Aborting")
            exit(-1)

        params = {
            'lambda': self.lmbd,
            'g_A': self.g_A,
            'g_X': self.g_X
        }
        return params

    @staticmethod
    def combinator(g1, g2, lmbd):
        return g1 + g2 + lmbd * g1 * g2
