from distutils.core import setup
import py2exe

setup(console=['generate.py'],
      options={"py2exe":{"includes":["pygments.styles.default"]}})