import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__),fname)).read()

setup(name='qeijo',
      version='0.1',
      description="Lightweight library to easily launch ab initio calculations with Quantum Espresso.",
      long_description=read('Readme.MD'),
      requires=['subprocess','io','os','shlex','exceptions'],
      author='Roberto Gomes de Aguiar Veiga',
      url="https://github.com/rgaveiga/qeijo",
      packages=['qeijo'])


