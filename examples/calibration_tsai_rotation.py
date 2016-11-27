# create the equivalent to tsai 1999 problem
# given the four frames: b,e,c,m
# we look for the fixed transform ec, while the others do change
#
# this becomes equivalent of solving following: Ai X = X Bi
# so we generate X, A and then obtain B as in: inv(X) A X = B
# feed solver for multiple pairs of (A,B)
#
#
# NOTE: this can be done also in SE(3) nothing that the two problems are separate, first SO(3) and then
# work on the euclidean part
import autograd.numpy as np

from pymanopt import Problem
from pymanopt.solvers import TrustRegions
from pymanopt.manifolds import Rotation

if not hasattr(__builtins__, "xrange"):
    xrange = range

if __name__ == "__main__":

	# generate X
	# generate Ai
	# compute Bi
	# solve problem
	k = 10
	manifold = Rotation(3,1)
	manifoldk = Rotation(3,k) # for simpliyfing A generation

	X = manifold.rand()
	iX = manifold.inverse(X)
	A = manifoldk.rand()
	B = np.zeros(np.shape(A)) # TODO same shape of A
    for i in xrange(k):
    	Ai = A[i]
    	B[i] = np.dot(iX,Ai,X)

    def cost(Xp):
        # AX = XB measured with the proper distance
        l = np.dot(A,Xp)
        r = np.dot(Xp,B)
        return manifoldk.dist(l,r)

    # first-order, second-order
    solver = TrustRegions()

    problem = Problem(manifold=manifold, cost=cost, verbosity=0)
    Xsol = solver.solve(problem)

    print("sol\n",Xsol)
    print("ref\n",X)
