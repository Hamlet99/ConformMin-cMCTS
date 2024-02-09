import numpy as np
from random import random


class Perturbation:
    def __init__(self, max_mutation):
        self.max_mutation = max_mutation

    @staticmethod
    def scale(depth, a, max_depth):
        """
        Scale the mutation based on the depth of the tree.
        :param depth:  depth of the tree.
        :type depth:  int
        :param a:  scaling parameter
        :type a:  float
        :param max_depth:  maximum depth of the tree
        :type max_depth:  int
        :return: scaling factor
        :rtype:  float
        """
        if random() > 0.2:
            pass
        else:
            a = 0

        depth_scale = np.exp(-a * (depth / max_depth) ** 2)
        return depth_scale

    def perturb(self, parameters, depth, a, max_depth):
        """
        Perturb the parameters by a random amount.
        :param parameters:  list of parameters to be perturbed
        :type parameters:  list
        :param depth:  depth of the tree
        :type depth:  int
        :param a:  scaling parameter
        :type a:  float
        :param max_depth:  maximum depth of the tree
        :type max_depth:  int
        :return:  list of perturbed parameters
        :rtype:  list
        """

        x = parameters.copy()
        u = np.random.normal(0.0, 360.0, len(x))
        delta = (
                (u / np.linalg.norm(u)) * self.max_mutation * self.scale(depth, a, max_depth)
        )
        x += delta
        x[x > 360] = 360
        x[x < 0] = 0

        return x
