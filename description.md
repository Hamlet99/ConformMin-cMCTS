## Monte Carlo Tree Search in Continuous Action Space

Monte Carlo Tree Search (MCTS) is a widely employed global optimization algorithm, particularly well-known for its effectiveness in discrete search problems, prominently in gaming domains such as AlphaGo. Originating from the fusion of tree search algorithms with principles of reinforcement learning, MCTS offers a heuristic and probabilistic approach to decision-making. Initially conceived for board games, MCTS constructs a decision tree reflecting potential moves, executing a sequence of fundamental steps:

1. **Selection**: Navigate down the tree to a leaf node.
2. **Expansion**: Add child nodes for unexplored moves.
3. **Simulation**: Conduct random playout from the new node.
4. **Backpropagation**: Update node statistics based on playout outcome.

In Continuous Monte Carlo Tree Search (c-MCTS), adapted for non-discrete problems without a clear "win" condition operating within non-discrete spaces, each node represents a point in the global search space. Nodes act as anchor points, and the decision tree refines the total search space. 

Exploration and exploitation are crucial strategies in c-MCTS, where the algorithm simultaneously explores potential pathways and exploits pathways with the highest estimated value.

### Enhanced Exploration and Avoiding Degeneracy

Continuous action space often faces the challenge of degeneracy, where different branches converge to the same area, leading to inefficient exploration. To address this, a uniqueness criteria function is introduced, scaling down exploration for degenerate solutions and scaling up for under-explored regions:

$$\ f\( \vec r_i\) = {1.5 \over 1 + \sum_{j \neq i}^{N_{points}} \delta \(\left\lvert r_{i} - r_{j}\right\rvert\)}$$

$$\ \delta \(\left\lvert r_{i} - r_{j}\right\rvert\) = \begin{cases} 1, &  \left\lvert r_i - r_j\right\rvert\ < r_{max}, \\
0, & \left\lvert r_i - r_j\right\rvert\ \ge r_{max}  \end{cases} $$


where $\vec{r_i}$ is the position of the $i_{th}$ point, $N_{points}$ is the total number of points, and $r_{max}$ is the maximum distance between points. The function scales exploration down towards zero for degenerate solutions and up for unique or under-explored regions, promoting efficient search across diverse solution spaces.

The node selection rule, Upper Confidence Bound for Parameters (UCP), integrates this function to guide exploration effectively:


$$\ UCP(\theta_j) = -min(r_1, r_2, ..., r_{n_i}) + c \cdot f(\theta_j) \cdot \sqrt { lnN_i\over n_i }$$

where $\theta_j$ represents node $j$ in the c-MCTS structure, $r$ is the reward for a given playout, $c$ is the exploration constant, $f(\theta_j)$ is the uniqueness criteria value for this node, $i$ is the number of playout samples taken by this node and all of its child nodes, and $N_i$$
is a similar value as $n_i$ except it is the parent node’s playout count instead of this node’s.

### Playouts with Adaptive Sampling

Playouts involve random sampling around a given point rather than clear "win" conditions in continuous action space. Playouts are executed by performing random vector displacements from the parameter set contained in the node, guided by the c-MCTS algorithm.  Attention is paid to the stochastic sampling process to ensure the generation of highly-representative sample points, overcoming challenges posed by the curse of dimensionality.

### Exploitation in Continuous Action Space

Efficient convergence in continuous space requires controlling the step size to avoid randomness or overshooting. A window scaling scheme is used, adjusting the maximum vector length based on the depth of the node in the c-MCTS tree. This incremental refinement of the phase space restores correlation between parent and child nodes, ensuring logical exploration and faster convergence.

For more details, refer to the 3rd section of the [supplementary materials](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27849-6/MediaObjects/41467_2021_27849_MOESM1_ESM.pdf) **Manna S. et al**.

