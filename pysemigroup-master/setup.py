from setuptools import setup
setup(name='pysemigroup',
      version='0.3b4',
      description='A tool to manipulate transitions semigroups and display them',
      url='https://gitlab.inria.fr/cpaperma/pysemigroup',
      author='Charles Paperman',
      author_email='charles.paperman@univ-lille.fr',
      license='GPLv2+',
      packages=['pysemigroup'],
      keywords=' automata semigroups monoid regular language',
      install_requires=["networkx","numpy"])
