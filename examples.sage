from F5 import *

def example(n, delta):
    macaulay2.eval("""R := GF(9001)[vars(0.."""+str(n)+""")];""");
    macaulay2.eval("""F := {};""");
    macaulay2.eval("""for i from 1 to """+str(n)+"""  do F = append(F, random("""+str(delta)+""",R))""");
    R = macaulay2('R').sage()
    F = macaulay2('F').sage()
    G_F5 = F5(F,n*delta-n+1)
    gb = list(G_F5[0][-1].values())
    print("Polynomial system: ", F[1:])
    print("Groebner basis: ", gb)
    print("Number of operations: ", G_F5[-2])
    print("Total polynomials added: ", len(gb))
    print("Basis is Groebner? ", R.ideal(gb).basis_is_groebner())
