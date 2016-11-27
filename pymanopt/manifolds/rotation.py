import numpy as np
from scipy.linalg import expm,logm

from pymanopt.tools.multi import multiprod, multitransp, multisym, multiskew
from pymanopt.manifolds.manifold import Manifold

if not hasattr(__builtins__, "xrange"):
    xrange = range


class Rotation(Manifold):
    """
        Rotation or Product of rotation SO(n)
    """

    def __init__(self, n, k=1):
        # Check that n is greater than or equal to p
        if n < 1:
            raise ValueError("Need n >= 1. Values supplied were n = %d " % n)
        if k < 1:
            raise ValueError("Need k >= 1. Value supplied was k = %d." % k)
        # Set the dimensions of the Rotation matrix
        self._n = height
        self._k = k
        self._dim = self._k * self._n * self._n

    @property
    def dim(self):
        return self._dim

    def __str__(self):
        if self._k == 1:
            return "Rotation manifold SO(%d)" % (self._n)
        elif self._k >= 2:
            return "Product Rotation manifold SO(%d)^%d" % (
                self._n, self._k)

    @property
    def typicaldist(self):
        return np.sqrt(self._n * self._k)

    def dist(self, X, Y):
        # Geodesic distance on the manifold
        raise np.linalg.norm(self.log(X, Y))

    def inverse(self,X);
        return X.T

    def inner(self, X, G, H):
        # Inner product (Riemannian metric) on the tangent space
        # For the Rotation this is the Frobenius inner product.
        return np.tensordot(G, H, axes=G.ndim)

    def proj(self, X, U):
        return multiskew(multiprod(multitransp(X), U))

    egrad2rgrad = proj

    # if det < 0 then flip some columns
    @staticmethod
    def flipqrsign(q,r):
        return np.dot(q, np.diag(np.sign(np.sign(np.diag(r))+.5)))

    def ehess2rhess(self, X, egrad, ehess, H):
        # Convert Euclidean into Riemannian Hessian.
        XtG = multiprod(multitransp(X), egrad)
        symXtG = multisym(XtG)
        XtEhess = multiprod(Xt, symXtG)
        return multiskew( XtEhess - multiprod(H, symXtG) )

    # (a retraction Rx is a mapping from the tangent space TxM to the manifold M at x onto M).
    def retr(self, X, U):
        XNew = X + multiprod(X, U);
        if self._k == 1:
            # Calculate 'thin' qr decomposition of X + G
            q, r = np.linalg.qr(Xnew)
            # Unflip any flipped signs
            return Rotation.flipqrsign(q,r)
        else:
            for i in xrange(self._k):
                q, r = np.linalg.qr(XNew[i])
                XNew[i] = Rotation.flipqrsign(q,r)
        return XNew

    def retr2(self, X, U):
        XNew = X + multiprod(X, U);
        if self._k == 1:
            Uk, ignore, Vk = np.linalg.svd(XNew)
            XNew = np.dot(Uk,Vk)
        else:
            for i in xrange(self._k):
                Uk, ignore, Vk = np.linalg.svd(XNew[i])
                XNew[i] = np.dot(Uk,Vk)
        return XNew


    def norm(self, X, G):
        # Norm on the tangent space of the Rotation is simply the Euclidean
        # norm.
        return np.linalg.norm(G)

    # Generate random Rotation point using qr of random normally distributed
    # matrix.
    # If Q is in O(n) but not in SO(n), we permute the two first
    # columns of Q such that det(new Q) = -det(Q), hence the new Q will
    # be in SO(n), uniformly distributed.
    def rand(self):
        if self._k == 1:
            A = np.random.randn(self._n, self._n)
            q, r = np.linalg.qr(A)
            return Rotation.flipqrsign(q,r)
        else:
            # Generated as such, Q is uniformly distributed over O(n), the set
            # of orthogonal matrices.
            X = np.zeros((self._k, self._n, self._n))
            for i in xrange(self._k):
                # replace with randn n n k
                A = np.random.randn(self._n, self._n)
                q, r = np.linalg.qr(A)
                X[i] = Rotation.flipqrsign(q,r)
            return X    

    @staticmethod
    def randskew(n, k):
        # TODO
        S = np.zeros((k,n,n))
        # only for L, others are 0
        S(L) = randn(size(L));      
        S = S - multitransp(S);  
        #[I J] = find(triu(ones(n), 1));
        #K = repmat(1:N, n*(n-1)/2, 1);
        #L = sub2ind([n n N], repmat(I, N, 1), repmat(J, N, 1), K(:));

    def randvec(self, X):
        U = Rotation.randskew(n, k);
        U = U / np.linalg.norm(U)
        return U

    def transp(self, x1, x2, d):
        return d

    # performs logarithm of inv(X) Y and return it in matrix form
    def log(self, X, Y):
        U = multiprod(multitransp(X), Y);
        # TODO multilog, but multilog is errimpl
        for i in xrange(self._k):
            U[i] = real(logm(U[i]));
        return multiskew(U) 

    def exp(self, X, U):
        if self._k == 1:
            Y = np.dot(X,expm(U))
        else:
            exptU = np.zeros(np.shape(X))
            # TODO multiexp, but multiexp is errimpl
            for i in xrange(self._k):
                exptU[i] = expm(U[i])
            Y = multiprod(X, exptU)
        return Y

    def pairmean(self, X, Y):
        raise NotImplementedError
