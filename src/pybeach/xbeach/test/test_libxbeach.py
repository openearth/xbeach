import os
import sys
import copy
import collections
try:
    import unittest2 as unittest
except ImportError:
    import unittest
from nose import with_setup
from .. import libxbeach
import numpy as np
import logging
#logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel( logging.DEBUG )
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s - %(name)s: %(message)s'))
logger.addHandler(handler)



import matplotlib.pyplot as plt
plt.interactive(True)
import time

dllsuffix = collections.defaultdict(lambda:'.so')
dllsuffix['darwin'] = '.dylib'
dllsuffix['win32'] = '.dll'
dllsuffix['win64'] = '.dll'
XBEACHLIB = os.path.join(
    os.path.dirname(__file__),
    '../../../xbeachlibrary/.libs/libxbeach.0' + dllsuffix[sys.platform]
    )
WD = os.path.join(
    os.path.dirname(__file__),
    '../../../../../branches/rewind/data/example1'
    )
class TestXBeach(unittest.TestCase):
    def setUp(self):
        "set up test fixtures"
        self.libpath = XBEACHLIB
        self.workingdir = WD
        self.xb = libxbeach.XBeach(
            libpath=self.libpath,
            workingdir=self.workingdir
            )
        self.xb.init()
    def tearDown(self):
        "tear down test fixtures"
        self.xb.finalize()
    def test_get_nparameter(self):
        self.assertGreaterEqual(self.xb.get_nparameter(), 224)
    def test_get_parameternamebyindex(self):
        self.assertEqual(self.xb.get_parameternamebyindex(0), 'depfile')
        self.assertEqual(self.xb.get_parameternamebyindex(20), 'disch_timeseries_file')
    def test_get_parameternames(self):
        names = self.xb.get_parameternames()
        self.assertGreaterEqual(len(names), 224)
        self.assertTrue('disch_loc_file' in names)
    def test_get_parametertype(self):
        name = 'eps'
        self.assertEqual(self.xb.get_parametertype(name), 'r')
        name = 'morfacopt'
        self.assertEqual(self.xb.get_parametertype(name), 'i')
        name = 'depfile'
        self.assertEqual(self.xb.get_parametertype(name), 'c')
    def test_get_parameter(self):
        name = 'eps'
        self.assertEqual(self.xb.get_parameter(name), 0.005)
        name = 'morfacopt'
        self.assertEqual(self.xb.get_parameter(name), 1)
        name = 'depfile'
        self.assertEqual(self.xb.get_parameter(name), 'z.grd')
    def test_get_parameters(self):
        parameters = self.xb.get_parameters()
        self.assertGreaterEqual(len(parameters), 224)
    def test_set_parameter(self):
        # Setting parameters should not work unless the parameters are not actual paramters (t, tnext)
        self.assertEqual(self.xb.get_parameter('g'), 9.81)
        self.assertRaises(ValueError, self.xb.set_parameter, 'g', 10.0)
        self.assertEqual(self.xb.get_parameter('g'), 9.81)
        
        self.assertEqual(self.xb.get_parameter('nx'), 707)
        self.assertRaises(ValueError, self.xb.set_parameter, 'nx', 1)
        self.assertEqual(self.xb.get_parameter('nx'), 707)

        self.assertGreaterEqual(self.xb.get_parameter('t'), 0.0)
        self.assertEqual(self.xb.set_parameter('t', 1.0), None)
        self.assertEqual(self.xb.get_parameter('t'), 1.0)
        self.assertEqual(self.xb.set_parameter('t', 0.0), None)
    def test_get_arraynames(self):
        self.assertGreaterEqual(len(self.xb.get_arraynames()), 30)
    def test_get_arraytype(self):
        self.assertEqual(self.xb.get_arraytype('zb'), 'r')
        self.assertEqual(self.xb.get_arraytype('wetz'), 'i')
    def test_get_array(self):
        self.assertEqual(self.xb.get_array('zb')[0,0], -19.96)
        self.assertEqual(self.xb.get_array('wetz')[0,0], 1)
        self.assertEqual(self.xb.get_array('dy'), 0.0)
        self.assertEqual(self.xb.get_array('dx'), 0.0)
    def test_get_arraydimsize(self):
        self.assertEqual(self.xb.get_arraydimsize('zb',0), 708)
        self.assertEqual(self.xb.get_arraydimsize('zb',1), 1)
    def test_get_arraydimshape(self):
        self.assertEqual(self.xb.get_arrayshape('zb'), (708,1))
    def test_get_arrays(self):
        arrays = self.xb.get_arrays()
        self.assertGreaterEqual(len(arrays), 30)
    def test_set_arrays(self):
        arrays = self.xb.get_arrays()
        for name in arrays:
            self.xb.set_array(name, arrays[name])
    def test_executestep(self):
        self.xb.executestep()
        self.assertGreater(self.xb.get_parameter('t'), 0.0)


class IntegrationTest(unittest.TestCase):
    def setUp(self):
        "set up test fixtures"
        self.libpath = XBEACHLIB
        self.workingdir = WD
        self.xb = libxbeach.XBeach(
            libpath=self.libpath,
            workingdir=self.workingdir
            )
        self.xb.init()
    def tearDown(self):
        "tear down test fixtures"
        self.xb.finalize()
    def test_rewind(self):
        # arrays point directly to XBeach memory....
        # Do a first timestep
        logger.debug('step')
        self.xb.executestep()
        t0_arrays = copy.deepcopy(self.xb.get_arrays())
        told = self.xb.get_parameter('t')
        logger.debug('t0 = %s', self.xb.get_parameter('t'))

                
        tnext = 10
        self.xb.set_parameter('tnext', tnext)
        while self.xb.get_parameter('t') < tnext:
            self.xb.executestep()
            self.xb.output()
            self.xb.set_parameter('tnext', tnext)
        t1a = self.xb.get_parameter('t')
        logger.debug('t1 = %s', self.xb.get_parameter('t'))
        t1a_arrays = copy.deepcopy(self.xb.get_arrays())
        
        self.xb.set_parameter('t', told)
        for name in t0_arrays:
            self.xb.set_array(name, t0_arrays[name])
        logger.debug('t2 = %s', self.xb.get_parameter('t'))
            
        self.xb.set_parameter('tnext', tnext)
        while self.xb.get_parameter('t') < tnext:
            self.xb.executestep()
            self.xb.output()
            self.xb.set_parameter('tnext', tnext)
        t1b = self.xb.get_parameter('t')
        logger.debug('t3 = %s', self.xb.get_parameter('t'))
        

        t1b_arrays = copy.deepcopy(self.xb.get_arrays())

        self.assertEqual(set(t1a_arrays), set(t1b_arrays))
        for name in t1a_arrays:
            if isinstance(t1a_arrays[name], np.ndarray) and t1a_arrays[name].size > 0:
                maxdiff = np.abs(t1a_arrays[name]- t1b_arrays[name]).max()
            else:
                maxdiff = np.abs(t1a_arrays[name]- t1b_arrays[name])
            if maxdiff > 0.00001:
                print name, maxdiff
        
