import numpy as np
from math import log


class Node:
    """
    Represents a node in the Monte Carlo Tree Search (MCTS) algorithm. Each node contains information about its
    index, depth, parent, number of visits: the number of times the node has been visited, children: list of indices
    of child nodes, and playout data: list of indices representing playout data associated with the node.

    :param idx:  The index of the node.
    :type idx:  int
    :param depth:  The depth of the node in the search tree.
    :type depth:  int
    :param parent:  The index of the parent node. Default is -1.
    :type parent:  int or None
    """

    def __init__(self, idx, depth, parent=-1):
        """
        Constructor method
        """
        self.idx = idx
        self.depth = depth
        self.parent = parent
        self.visits = 1
        self.children = []
        self.playout_data = []


class cMCTS:
    """
    Implements the Monte Carlo Tree Search (MCTS) algorithm.

    :param root_data:  The initial randomly generated parameter list for the root node.
    :type root_data:  list
    :param perturbate:  The perturbation function to perturb the data. Output is a list of perturbed parameters.
    :type perturbate:  function
    :param evaluate:  The evaluation function to evaluate the data. Output is a list: [list of parameters, score].
    :type evaluate:  function
    :param n_iterations:  The number of iterations for the MCTS algorithm. Default is 200.
    :type n_iterations:  int
    :param n_expand:  The number of expansions to perform during each iteration. Default is 10.
    :type n_expand:  int
    :param n_playouts:  The number of playouts to perform. Default is 20.
    :type n_playouts:  int
    :param explore_constant:  The exploration constant for the Upper Confidence Bound (UCB) formula. Default is 1.0.
    :type explore_constant:  float
    :param max_depth:  The maximum depth of the search tree. Default is 12.
    :type max_depth:  int
    :param a:  A scaling parameter used for perturbation. Default is 3.0.
    :type a:  float
    :param selected_node:  The currently selected node in the MCTS algorithm. Default is 0.
    :type selected_node:  int
    :param nodes:  A dictionary containing Node objects with their respective indices as keys.
    :type nodes:  dict
    :param score_list:  A list containing scores of evaluated data. Default is [].
    :type score_list:  list
    :param parameter_list:  A list containing all evaluated parameters. Default is [].
    :type parameter_list:  list
    """

    def __init__(self, root_data, perturbate, evaluate, **kwargs):
        """
        Constructor method
        """
        self.root_data = root_data
        self.perturbate = perturbate
        self.evaluate = evaluate
        self.n_iterations = kwargs.get('n_iterations', 200)
        self.n_expand = kwargs.get('n_expand', 10)
        self.n_playouts = kwargs.get('n_playouts', 20)
        self.explore_constant = kwargs.get('explore_constant', 1)
        self.max_depth = kwargs.get('max_depth', 12)
        self.a = kwargs.get('a', 3)

        # Initialize the root node
        self.selected_node = 0
        self.nodes = {0: Node(0, 0, None)}
        root_params, root_energy = evaluate(self.root_data)  # must be adapted for molecule torsions
        self.score_list = [root_energy]
        self.parameter_list = [root_params]

    def playouts(self, idx):
        """
        Perform playouts for a given node. Perturb the data, evaluate it, and add the data and
        score to the respective lists.

        :param idx:  The index of the node.
        :type idx:  int
        :return:  The updated index.
        :rtype:  int
        """

        initial_node = idx
        for i in range(self.n_playouts):
            idx += 1
            playdata = self.perturbate(self.parameter_list[initial_node], depth=self.nodes[initial_node].depth,
                                       a=self.a, maxdepth=self.max_depth)
            playdata_params, play_energy = self.evaluate(playdata)
            self.score_list.append(play_energy)
            self.parameter_list.append(playdata_params)
            self.nodes[initial_node].playout_data.append(idx)

        return idx

    def expansion_simulation(self, parent_node, idx):
        """
        Perform the expansion and simulation steps of the MCTS algorithm. Perturb the data, evaluate it, and
        add the data and energy to the respective lists. Update the index and the selected node.

        :param parent_node:  The index of the parent node.
        :type parent_node:  int
        :param idx:  The index of the node.
        :type idx:  int
        :return: The updated index.
        :rtype:  int
        """
        play_indexes = self.nodes[parent_node].playout_data
        play_energies = [self.score_list[i] for i in play_indexes]
        best_play_index = play_indexes[play_energies.index(min(play_energies))]

        idx += 1
        self.nodes[parent_node].visits += 1
        self.nodes[parent_node].children.append(idx)

        depth = self.nodes[parent_node].depth + 1
        self.nodes.update({idx: Node(idx, depth, parent_node)})

        data = self.perturbate(self.parameter_list[best_play_index], depth=depth, a=self.a, maxdepth=self.max_depth)
        parameters, energy = self.evaluate(data)
        self.score_list.append(energy)
        # depending on the application, parameters might be the same as perturbed data
        self.parameter_list.append(parameters)

        idx = self.playouts(idx)

        return idx

    def get_branch(self, node):
        """
        Get the branch of a given node.

        :param node:  The index of the node.
        :type node:  int
        :return:  List of indices representing the branch from the root node to the given node.
        :rtype:  list
        """
        list_branch_indexes = []
        if node in self.nodes.keys():
            list_branch_indexes += self.nodes[node].children
            for child in self.nodes[node].children:
                list_branch_indexes += self.get_branch(child)
        return list_branch_indexes

    def backpropagation_selection(self):
        """
        Back-propagate and select the best node based on the Upper Confidence Bound (UCB) formula.

        :return: The index of the best node.
        :rtype:  int
        """

        selection_energies = []
        nodes = [i for i, j in self.nodes.items() if j.parent != -1]

        for n in nodes:
            parent = self.nodes[n].parent
            depth = self.nodes[n].depth
            if parent is None or depth > self.max_depth:
                selection_energies.append(-1e300)

            else:
                branch = self.get_branch(n)
                all_indexes = []
                for i in [n] + branch:
                    all_indexes += self.nodes[i].playout_data

                best_reward = min([self.score_list[i] for i in branch + all_indexes])

                UCB_score = -best_reward + self.explore_constant * np.sqrt(
                    log(self.nodes[parent].visits) / self.nodes[n].visits)

                selection_energies.append(UCB_score)

        c = zip(selection_energies, nodes)
        c = sorted(c, key=lambda x: -x[0])
        selection_energies, nodes = zip(*c)

        return nodes[0]

    def run(self):
        """
        Run the MCTS algorithm. Perform playouts, expansion, simulation, and back-propagation for the selected node and
        prints the total number of evaluations and the current best score.

        :return: None
        """

        idx = self.playouts(0)  # run playouts for the first node

        for _ in range(self.n_iterations):
            for _ in range(self.n_expand):
                idx = self.expansion_simulation(self.selected_node, idx)
                self.selected_node = self.backpropagation_selection()

            print(f"Evaluations made: {len(self.score_list)}, Best energy yet: {min(self.score_list)}")
