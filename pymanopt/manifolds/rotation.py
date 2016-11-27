import numpy as np
from scipy.linalg import expm

from pymanopt.tools.multi import multiprod, multitransp, multisym, multiskew
from pymanopt.manifolds.manifold import Manifold

if not hasattr(__builtins__, "xrange"):
    xrange = range


class Rotation(Manifold):
    """
    """

    def __init__(self, n, k=1):
        # Check that n is greater than or equal to p
        if n < 1:
            raise ValueError("Need n >= 1. Values supplied were n = %d " % n)
        if k < 1:
            raise ValueError("Need k >= 1. Value supplied was k = %d." % k)
        # Set the dimensions of the Stiefel
        self._n = height
        self._k = k
        self._dim = self._k * self._n

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
        raise np.linalg.norm(X, self.log(X, Y))

    def inverse(self,X);
        return X.T

    def inner(self, X, G, H):
        # Inner product (Riemannian metric) on the tangent space
        # For the stiefel this is the Frobenius inner product.
        return np.tensordot(G, H, axes=G.ndim)

    def proj(self, X, U):
        return multiskew(multiprod(multitransp(X), U))

    egrad2rgrad = proj

    def ehess2rhess(self, X, egrad, ehess, H):
        # Convert Euclidean into Riemannian Hessian.
        XtG = multiprod(multitransp(X), egrad)
        symXtG = multisym(XtG)
        XtEhess = multiprod(Xt, symXtG)
        return multiskew( XtEhess - multiprod(H, symXtG) )

    # SAME
    def retr(self, X, G):
        XNew = X + multiprod(X, tU);
        if self._k == 1:
            # Calculate 'thin' qr decomposition of X + G
            q, r = np.linalg.qr(X + G)
            # Unflip any flipped signs
            XNew = np.dot(q, np.diag(np.sign(np.sign(np.diag(r))+.5)))
        else:
            XNew = X + G
            for i in xrange(self._k):
                q, r = np.linalg.qr(XNew[i])
                XNew[i] = np.dot(q, np.diag(np.sign(np.sign(np.diag(r))+.5)))
        return XNew

    def retr2(self, X, G):
        XNew = X + multiprod(X, tU);
        if self._k == 1:
            Uk, ignore, Vk = np.linalg.svd(XNew)
            XNew = np.dot(Uk,Vk)
        else:
            for i in xrange(self._k):
                Uk, ignore, Vk = np.linalg.svd(XNew[:,k])
                XNew[:,:,k] = np.dot(Uk,Vk)
        return XNew


    def norm(self, X, G):
        # Norm on the tangent space of the Stiefel is simply the Euclidean
        # norm.
        return np.linalg.norm(G)

    # Generate random Stiefel point using qr of random normally distributed
    # matrix.
    def rand(self):
        # TODO randrot
        """
        R = zeros(n, n, N);
    
    for i = 1 : N
        
        % Generated as such, Q is uniformly distributed over O(n), the set
        % of orthogonal matrices.
        A = randn(n);
        [Q, RR] = qr(A);
        Q = Q * diag(sign(diag(RR))); %% Mezzadri 2007
        
        % If Q is in O(n) but not in SO(n), we permute the two first
        % columns of Q such that det(new Q) = -det(Q), hence the new Q will
        % be in SO(n), uniformly distributed.
        if det(Q) < 0
            Q(:, [1 2]) = Q(:, [2 1]);
        end
        
        R(:, :, i) = Q;
        
    end"""
        return
        if self._k == 1:
            X = np.random.randn(self._n, self._p)
            q, r = np.linalg.qr(X)
            return q

        X = np.zeros((self._k, self._n, self._p))
        for i in xrange(self._k):
            X[i], r = np.linalg.qr(np.random.randn(self._n, self._p))
        return X

    def randvec(self, X):
        U = randskew(n, k);
        U = U / np.linalg.norm(U)
        return U

    def transp(self, x1, x2, d):
        return d

    def log(self, X, Y):
        U = multiprod(multitransp(X), Y);
        for i in xrange(self._k):
            U[:,:,i] = real(logm(U[:, :, i]));
        return multiskew(U) 

    def exp(self, X, U):
        if self._k == 1:
            Y = np.dot(X,expm(U))
        else:
            exptU = np.zeros(np.shape(X))
            for i in xrange(self._k):
                exptU[:,:,i] = expm(exptU[:,:,i])
            Y = multiprod(X, exptU)
        return Y

    def pairmean(self, X, Y):
        raise NotImplementedError
