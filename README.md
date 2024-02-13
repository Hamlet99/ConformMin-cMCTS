[![Documentation Status](https://readthedocs.org/projects/conformmin-cmcts/badge/?version=latest)](https://conformmin-cmcts.readthedocs.io/en/latest/?badge=latest)
## ConformMin-cMCTS
**ConformMin-cMCTS**: Repository for conformational energy global minimization using continuous Monte Carlo Tree Search (cMCTS) algorithm.

Continuous Monte Carlo Tree Search (cMCTS) algorithm was implemented based on papers listed below:
- Manna, S., Loeffler, T.D., Batra, R. et al. Learning in continuous action space for developing high dimensional potential energy models. Nat Commun 13, 368 (2022). https://doi.org/10.1038/s41467-021-27849-6
- Banik, S et al. A Continuous Action Space Tree search for INverse desiGn (CASTING) framework for materials discovery. npj Comput Mater 9, 177 (2023). https://doi.org/10.1038/s41524-023-01128-y


## Installation
To install the package, you can use `pip` after cloning the repository to your local machine. 
```bash
git clone git@github.com:Hamlet99/ConformMin-cMCTS.git
cd /path/to/ConformMin-cMCTS/
pip install -e .
```

## Usage

### Specify your Model

Specify all parameters of the model to be solved in a `config_mcts.json` configuration file, with the following structure:

```json
{
    "n_iterations": 200,
    "n_expand": 10,
    "n_playouts": 20,
    "explore_constant": 1.0,
    "max_depth": 12,
    "a": 3
}
```
- `"n_iterations"` (int) The number of iterations for the MCTS algorithm. Default is 200.
- `"n_expand"` (int)  The number of expansions to perform during each iteration. Default is 10.
- `"n_playouts"` (int) The number of playouts to perform. Default is 20.
- `"explore_constant"` (float) The exploration constant for the Upper Confidence Bound (UCB) formula. Default is 1.0.
- `"max_depth"` (int) The maximum depth of the search tree. Default is 12.
- `"a"` (float) A scaling parameter used for perturbation. Default is 3.0.

### Run the MCTS algorithm

To run the MCTS algorithm, you can use the following code in your Python script or Jupyter notebook (`example.ipynb`):


``` python
import environment
from cMCTS import mcts
from cMCTS import perturbation

smiles = 'c1cnccc1c1ccncc1NC'

#Initialize the environment and draw the marked molecule using RDKit
env = environment.Environment(smiles)
root_data = env.create_root_parameters()
env.draw()

#Initialize the perturbator
max_mutation = 10 #maximum mutation allowed
perturbator = perturbation.Perturbation(max_mutation)

# read the config file
with open('config_mcts.json') as f:
    config = json.load(f)
    print(config)

# set perturb and evaluate functions 
perturb = perturbator.perturb
evaluate = env.evaluate

# initialize, run cMCTS, and generate "mcts_results.json" file
optimizer = mcts.cMCTS(root_data, perturb, evaluate, **config)
optimizer.run()
optimizer.write()

# to visualize the optimization process use the following command  
optimizer.plot_minimization(lim_iterations=[0,1000])
```

License
ConformMin-cMCTS is released under the MIT License. See `LICENSE` for details.
