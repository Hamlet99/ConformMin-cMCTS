from setuptools import setup

setup(name='ConformMin-cMCTS',
      version='0.1',
      description="ConformMin-cMCTS: Repository for conformational energy minimization using continuous Monte Carlo Tree Search (cMCTS) algorithm.",
      author="Hamlet Khachatryan",
      author_email="hamlet.khachatryan@dtc.ox.ac.uk",
      packages=['cMCTS'],
      install_requires=["numpy", "pandas", "matplotlib", "sphinx", "pytest", "rdkit", "setuptools"],
      license='MIT',
      entry_points={
            "console_scripts": ["checkpoint = checkpoints:main"]
      }
)