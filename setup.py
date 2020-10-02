from setuptools import setup

setup(
    name='VETA',
    version='0.1',
    licence='GPL-3.0',
    author='Pedro Barbosa',
    author_email='psbpedrobarbosa@gmail.com',
    packages=['src'],
    entry_points={
        'console_scripts': ['veta=src.veta:main']
    }
)

