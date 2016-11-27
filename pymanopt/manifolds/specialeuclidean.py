from pymanopt.manifolds.manifold import Manifold,Product,Rotation,Euclidean

if not hasattr(__builtins__, "xrange"):
    xrange = range


class SpecialEuclidean(Product):

    def __init__(self, n, k=1):
    	# TODO: (R,e) (R,e) ... and not (R,...R) (e,...,e)
    	Product.__init__([Rotation(n,k),Euclidean(n,k)])