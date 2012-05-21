import os
import sys
import collections
from numbers import Number
from ctypes import c_void_p, c_char_p, c_int, c_double, c_char
from ctypes import POINTER, pointer
from ctypes import CDLL, string_at, addressof, byref,create_string_buffer
from numpy.ctypeslib import ndpointer, as_array
from numpy import float64, zeros, array, int32, ndarray
import subprocess

import logging
#logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel( logging.DEBUG )
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s - %(name)s: %(message)s'))
logger.addHandler(handler)


dllsuffix = collections.defaultdict(lambda:'.so')
dllsuffix['darwin'] = '.dylib'
dllsuffix['win32'] = '.dll'
dllsuffix['win64'] = '.dll'

def dlclose(handle):
    name = 'libdl' + dllsuffix[sys.platform]
    libdl = CDLL(name)
    libdl.dlerror.restype = c_char_p
    libdl.dlclose.argtypes = [c_void_p]
    logger.debug('Closing dll (%x)',handle)
    rc = libdl.dlclose(handle)
    if rc!=0:
        logger.debug('Closing failed, looking up error message')
        error = libdl.dlerror()
        logger.debug('Closing dll returned %s (%s)', rc, error)
        if error == 'invalid handle passed to dlclose()':
            raise ValueError(error)
    else:
        logger.debug('Closed')


def isloaded(lib):
    """return true if library is loaded"""
    libp = os.path.abspath(lib)
    ret = os.system("lsof -p %d | grep %s > /dev/null" % (os.getpid(), libp))
    return (ret == 0)
# We might want to run this part in a proxy....
class XBeach:
    """Proxy to the libxbeach library"""
    def __init__(self, libpath, workingdir=None):
        self.olddir = os.getcwd()
        if workingdir is None:
            self.workingdir = os.getcwd()
        else:
            self.workingdir = workingdir
        os.chdir(self.workingdir)
        self.libpath = os.path.abspath(libpath)
        self._lib = CDLL(self.libpath)
        self.shouldinitialize = True
    def init(self):
        """Initialise XBeach (only run once)"""
        logger.debug('Initializing...')
        if self.shouldinitialize:
            self._lib.init()
        self.shouldinitialize = False
    def executestep(self):
        """Execute a timestep, stops before or on tnext"""
        self._lib.executestep()
    def output(self):
        """Call XBeach output"""
        self._lib.outputext()
    def get_nparameter(self):
        n = c_int()
        self._lib.getnparameter(byref(n))
        return n.value
    def get_parameternamebyindex(self, i):
        # Just pass a pointer to get the address of the first character
        # All strings in XBeach are 1024 long...
        name = create_string_buffer(1024)
        index = c_int(i+1)
        length = c_int()
        self._lib.getparametername(byref(index), byref(name), byref(length))
        # Analogue to copybuffer...
        result = name.value
        return result
    def get_parameternames(self):
        parameternames = []
        for i in range(self.get_nparameter()):
            name = self.get_parameternamebyindex(i)
            parameternames.append(name)
        return parameternames
    def get_parameter(self, name):
        typecode = self.get_parametertype(name)
        c_name = create_string_buffer(name)
        namelength = c_int(len(name))
        valuelength = c_int()
        if typecode == 'c':
            value = c_char()
            self._lib.getcharparameter(byref(c_name), byref(value), byref(namelength), byref(valuelength))
            # extract the string at this location...
            result = string_at(addressof(value), valuelength)
        elif typecode == 'r':
            value = c_double()
            self._lib.getdoubleparameter(byref(c_name), byref(value), (namelength))
            result = value.value
        elif typecode == 'i':
            value = c_int()
            self._lib.getintparameter(byref(c_name), byref(value), (namelength))
            result = value.value
        else:
            raise ValueError('What is typecode {0} for variable {1}'.format(typecode, name))
        return result
    def get_parametertype(self, name):
        typecode = c_char()
        c_name = create_string_buffer(name)
        length = c_int(len(name))
        self._lib.getparametertype(byref(c_name), byref(typecode), byref(length))
        return typecode.value
    def get_parameters(self):
        parameters = {
            (name, self.get_parameter(name))
            for name
            in self.get_parameternames()
            }
        return parameters
    def set_parameter(self, name, value):
        typecode = self.get_parametertype(name)
        c_name = create_string_buffer(name)
        namelength = c_int(len(name))
        if typecode == 'r':
            value = c_double(value)
            code = self._lib.setdoubleparameter(byref(c_name), byref(value), (namelength))
            if (code != 0):
                raise ValueError('Error thrown by XBeach function: {0} for name {1}'.format('setdoubleparameter', c_name.value))
            result = value.value
        # if typecode == 'i':
        #     value = c_int(value)
        #     code = self._lib.setintparameter(byref(c_name), byref(value), (namelength))
        #     if (code != 0):
        #         raise ValueError('Error thrown by XBeach function: {0} for name {1}'.format('setintparameter', c_name.value))
        #     result = value.value
        else:
            raise ValueError('What is typecode {0} for variable {1}'.format(typecode, name))
    def get_narray(self):
        n = c_int()
        self._lib.getnarray(byref(n))
        return n.value
    def get_arraynamebyindex(self, i):
        # Just pass a pointer to get the address of the first character
        # All strings in XBeach are 1024 long...
        name = create_string_buffer(1024)
        index = c_int(i+1)
        length = c_int()
        self._lib.getarrayname(byref(index), byref(name), byref(length))
        # Analogue to copybuffer...
        result = name.value
        return result
    def get_arraynames(self):
        arraynames = []
        for i in range(self.get_narray()):
            name = self.get_arraynamebyindex(i)
            arraynames.append(name)
        return arraynames
    
    def get_arraytype(self, name):
        typecode = c_char()
        c_name = create_string_buffer(name)
        namelength = c_int(len(name))
        self._lib.getarraytype(byref(c_name), byref(typecode), byref(namelength))
        return typecode.value
    def get_arrayrank(self, name):
        rank = c_int()
        c_name = create_string_buffer(name)
        namelength = c_int(len(name))
        self._lib.getarrayrank(byref(c_name), byref(rank), byref(namelength))
        return rank.value
    def get_arraydimsize(self, name, dim):
        dimsize = c_int()
        dim = c_int(dim+1)
        c_name = create_string_buffer(name)
        namelength = c_int(len(name))
        self._lib.getarraydimsize(byref(c_name), byref(dim), byref(dimsize), byref(namelength))
        return dimsize.value
    def get_arrayshape(self, name):
        return tuple(
            self.get_arraydimsize(name, i)
            for i in range(self.get_arrayrank(name))
            )
            
            
    def get_array(self, name):
        typecode = self.get_arraytype(name)
        rank = self.get_arrayrank(name)
        c_name = create_string_buffer(name)
        namelength = c_int(len(name))
        result = None
        fun = None
        arraytype = None
        if rank == 0:
            shape = ()
        elif rank >= 1:
            shape = self.get_arrayshape(name)
        else:
            raise ValueError('Array {} not supported because {} not in (0, 2)'.format(name, rank))
        if typecode == 'c':
            raise ValueError('Getting character arrays is not supported')
        elif typecode == 'i':
            if rank > 2:
                raise ValueError('Getting arrays of rank {} is not supported'.format(rank))

            fun = getattr(self._lib, 'get{}dintarray'.format(rank))
            if rank == 0:
                arraytype = POINTER(c_int) # Double int pointer 
            else:
                arraytype = ndpointer(dtype=int32, shape=shape, ndim=rank, flags='F_CONTIGUOUS')
        elif typecode == 'r':
            if rank > 4:
                raise ValueError('Getting arrays of rank {} is not supported'.format(rank))

            fun = getattr(self._lib, 'get{}ddoublearray'.format(rank))
            if rank == 0:
                arraytype = POINTER(c_double)
            else:
                arraytype = ndpointer(dtype=float64, shape=shape, ndim=rank, flags='F_CONTIGUOUS')
        else:
            raise ValueError('Getting arrays of type {} is not supported'.format(typecode))
        fun.argtypes=[POINTER(c_name._type_), POINTER(arraytype), c_int]
        arrayp = arraytype()
        fun(c_name, byref(arrayp), namelength)
        if rank == 0:
            # Just return the number
            result = arrayp.contents.value
        else:
            result = array(arrayp)
        return result
    def get_arrays(self):
        arrays = {}
        for name in self.get_arraynames():
            arrays[name] = self.get_array(name)
        return arrays

    def set_array(self, name, value):
        typecode = self.get_arraytype(name)
        rank = self.get_arrayrank(name)
        c_name = create_string_buffer(name)
        namelength = c_int(len(name))
        fun = None
        arraytype = None
        
        shape = self.get_arrayshape(name)
        assert shape == array(value).shape, "Shapes not equal {} versus {} for variable {}".format(shape, value.shape, name)
        
        if typecode == 'c':
            raise ValueError('Setting character arrays is not supported')
        if rank > 4:
            raise ValueError('Getting arrays of rank {} is not supported'.format(rank))
        if typecode == 'r':
            dtype = float64 if rank > 0 else c_double
            typename = 'double'
        elif typecode == 'i':
            dtype = int32 if rank > 0 else c_double
            typename = 'int'
        fun = getattr(self._lib, 'set{}d{}array'.format(rank, typename))
        if rank == 0:
            arraytype = POINTER(dtype) 
        else:
            arraytype = ndpointer(dtype=dtype, shape=shape, ndim=rank, flags='F_CONTIGUOUS')
        fun.argtypes=[POINTER(c_name._type_), POINTER(arraytype), c_int]
        if isinstance(value, Number):
            arrayp = pointer(dtype(value))
        elif isinstance(value, ndarray):
            arrayp = value.ctypes.data_as(arraytype)
        else:
            raise ValueError('Expected value of type {} but got {}'.format((Number,ndarray), type(value)))
        fun(c_name, byref(arrayp), namelength)

    def finalize(self):
        logger.debug('Explicitly deleting library at {}'.format(self._lib))
        handle = self._lib._handle
        del self._lib
        logger.debug('Checking if %s is loaded %s', self.libpath, isloaded(self.libpath))
        while isloaded(self.libpath):
            dlclose(handle)
            logger.debug('Checking if %s is loaded %s', self.libpath, isloaded(self.libpath))
        # openfiles = subprocess.check_output(['lsof', '-p', str(os.getpid())])
        # logger.debug('Open files:\n%s',  openfiles)
        # reset the init method
        self.shouldinitialize = True
                             
