class SignedMatrix:
    # Matrix together with hashmap associating singature (index) to each row - has special rref function that respects signatures.
    def __init__(self, mat, sgn, d, parent):
        self.mat = mat
        self.signature = sgn
        self.d = d
        self.parent = parent

    # use position over term ordering
    def row_echelon_form_by_position(self):
        # returns a pair (M, n) where M is a new signed matrix which is the row-reduction of self via a sequence of
        # elementary row operations
        # keep track of number of operations
        num_operations = 0
        copy_mat = copy(self.mat)
        eliminated = True
        first_reduction = True
        # keep track of reductions
        rdxn = dict()
        for i in range(len(copy_mat.rows())):
            rdxn[i] = []
        while eliminated:
            eliminated = False
            for i, row in enumerate(copy_mat.rows()):
                for j in range(len(row)):
                    if row[j] != 0:
                        # j is the leading term of this row, so use it to kill everything with higher signature
                        for new_i, new_row in enumerate(copy_mat.rows()):
                            if new_row[j] != 0 and self.signature[i] < self.signature[new_i]:
                                # we can reduce
                                lam = -(new_row[j]/row[j])
                                copy_mat.add_multiple_of_row(new_i, i, lam)
                                eliminated = True
                                if first_reduction: # only count top-reductions
                                    num_operations += len(new_row)
                                rdxn[new_i].append((i,lam))
                        break
            first_reduction = False # stop counting arithmetic operations

        for i, row in enumerate(copy_mat.rows()):
            for j in range(len(row)):
                if row[j] != 0:
                    # j is the coefficient of the leading term of this row, so divide this row by it
                    copy_mat.rescale_row(i,1/row[j])
                    break

        return (SignedMatrix(copy_mat, self.signature, self.d, self.parent), num_operations, rdxn)

    # use term over position ordering
    def row_echelon_form_by_term(self):
        num_operations = 0
        copy_mat = copy(self.mat)
        eliminated = True
        first_reduction = True
        while eliminated:
            eliminated = False
            for i, row in enumerate(copy_mat.rows()):
                for j in range(len(row)):
                    if row[j] != 0:
                        # j is the leading term of this row, so use it to kill everything with higher signature
                        for new_i, new_row in enumerate(copy_mat.rows()):
                            if new_row[j] != 0 and self.signature[i][::-1] < self.signature[new_i][::-1]:
                                # we can reduce
                                lam = -(new_row[j]/row[j])
                                copy_mat.add_multiple_of_row(new_i, i, lam)
                                eliminated = True
                                if first_reduction:
                                    num_operations += len(new_row)
                        break
            first_reduction = False

        for i, row in enumerate(copy_mat.rows()):
            for j in range(len(row)):
                if row[j] != 0:
                    # j is the coefficient of the leading term of this row, so divide this row by it
                    copy_mat.rescale_row(i,1/row[j])
                    break

        return (SignedMatrix(copy_mat, self.signature, self.d, self.parent), num_operations)

    def add_row(self, f, index):
        # returns a new matrix which is self with a row added corresponding to polynomial f with signature index
        row = [f.monomial_coefficient(mon) for mon in self.monomials()]
        copy_mat = copy(self.mat)
        copy_signature = copy(self.signature)
        copy_mat = matrix(copy_mat.rows()+[row])
        copy_signature[copy_mat.nrows()-1] = index
        return SignedMatrix(copy_mat, copy_signature, self.d, self.parent)

    def monomials(self):
        # returns monomials of degree self.d in a list, sorted in decreasing order
        R = self.parent
        monomials = [R({tuple(a):1}) for a in WeightedIntegerVectors(self.d, [1]*R.ngens())]
        monomials.sort(reverse=True)
        return monomials

    def LT(self):
        # returns the leading terms of the (polynomials represented by) rows of self.mat
        monomials = self.monomials()
        leading_terms = []
        for row in self.mat.rows():
            for i in range(len(row)):
                if row[i] != 0:
                    leading_terms.append(monomials[i]*row[i])
                    break
        return set(leading_terms)

    def rows(self):
        # return set of (polynomials represented by) rows of self.mat
        monomials = self.monomials()
        r = []
        for row in self.mat.rows():
            polynomial = 0
            for j in range(len(row)):
                polynomial += row[j]*monomials[j]
            r.append(polynomial)
        return r 

def F5(F, D, order='position'):
    # F=(f_1,...,f_m) is a set of polynomials with degere d_1 <= d_2 <= ... <= d_m
    # D is maximal degree
    # returns the set of elements of degree at most D of reduced Grobner bases of (f_1,...,f_i) for each i
    operations = 0
    F.insert(0,0) # so that we can 1-index everything
    G = [{} for _ in range(len(F))] # initialize intermediate Grobner bases
    M = [[None for _ in range(len(F))] for _ in range(D+1)] # initialize Macaulay matrices
    M_red = [[None for _ in range(len(F))] for _ in range(D+1)] # initialize reduced Macaulay matrices
    sizes = [] # initialize list of sizes of Macaulay matrices
    rdxn = [[None for _ in range(len(F))] for _ in range(D+1)] # initialize list of reductions performed
    variables = list(F[1].parent().gens())
    variables.sort(reverse=True)
    for d in range(F[1].degree(),D+1):
        M[d][0] = SignedMatrix(matrix(QQ), dict(), d, F[1].parent())
        M_red[d][0] = SignedMatrix(matrix(QQ), dict(), d, F[1].parent())
        for i in range(1, len(F)):
            if d < F[i].degree(): M[d][i] = M[d][i-1] # Case 1: the degree of f_i is larger than d
            elif d == F[i].degree(): # Case 2: the degree of f_i is exactly d
                M[d][i] = M_red[d][i-1].add_row(F[i], (i,1))
            else: # Case 3: the degree of f_i is less than d
                M[d][i] = M_red[d][i-1]
                if M_red[d-F[i].degree()][i-1]:
                    Crit = M_red[d-F[i].degree()][i-1].LT() # build F_5 criterion list
                else:
                    Crit = []
                for row_num, sgn in [(r,s) for (r,s) in M[d-1][i].signature.items() if s not in M[d-1][i-1].signature.values()]:
                    _,u = sgn 
                    f = M[d-1][i].rows()[row_num]
                    if u == 1:
                        largest_var_in_u = 0
                    else:
                        largest_var_in_u = variables.index(u.variables()[-1]) # select which row to use to build new row
                    for j in range(largest_var_in_u,len(variables)):
                        if u*variables[j] not in Crit: # avoid signatures which F_5 criterion tells us are useless
                            M[d][i] = M[d][i].add_row(variables[j]*f, (i,u*variables[j]))

            # reduce Macaulay-like matrices
            if order == "position":
                M_red[d][i], op, rdxn[d][i] = M[d][i].row_echelon_form_by_position()
                operations += op
            else:
                M_red[d][i], op, rdxn[d][i] = M[d][i].row_echelon_form_by_term()
                operations += op
            sizes.append((M[d][i].mat.nrows(), M[d][i].mat.ncols()))
            # update Grobner bases
            for j,f in enumerate(M_red[d][i].rows()):
                if f.lt().reduce([g.lt() for g in G[i].values()]) != 0:
                    G[i][M_red[d][i].signature[j]] = f
    return (G,sizes,M[-1][-1].mat,M[-1][-1].signature, M_red, operations, rdxn)
