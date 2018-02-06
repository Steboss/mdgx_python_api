import numpy as np
from cpython cimport array

# Some notes
# - cython will auto convert simple C struct to python dict

cdef get_positions(MolecularDynamicsSystem mys):
    import numpy as np
    cdef int n_atoms = mys.crd.natom

    cdef double[:] arr = <double[:n_atoms*3]> mys.crd.loc

    return np.asarray(arr)

cdef set_positions(double[:] arr, MolecularDynamicsSystem mys):
    cdef int i
    for i in range(mys.crd.natom):
        mys.crd.loc[i] = arr[i]

cdef class setup:
    def __init__(self, prmtop, rst7):
        cdef PotentialFunction myu

        self._traj_control = create_trajcon_()
        self._myu = load_topology_(prmtop, &self._traj_control)
        self._mdsystem = create_mdsystem_(rst7, &self._myu)

        # trigger loading
        self.positions = self.positions

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.cleanup()

    def cleanup(self):
        destroy_mdsystem(&self._mdsystem)
        destroy_uform(&self._myu, &self._mdsystem)

    @property
    def n_atoms(self):
        return self._mdsystem.crd.natom

    def get_positions(self):
        # compat with pysander
        return self.positions

    def set_positions(self, arr):
        # compat with pysander
        self.positions = arr

    property positions:
        def __get__(self):
            return get_positions(self._mdsystem).reshape(self.n_atoms, 3)

        def __set__(self, arr):
            cdef double[:] arr_view = np.asarray(arr).flatten()
            load_coords_(&self._myu, &self._traj_control,
                          &arr_view[0], &self._mdsystem)
            InitExecon(&(self._mdsystem.etimers))
        
    def energy_forces(self):
        cdef Energy ene
        cdef cell* cell_ptr
        cdef int i, j, k
        cdef double[:] gradients = np.zeros(self.n_atoms*3)

        MMForceEnergy(&self._myu, &self._mdsystem, &self._traj_control)

        for i in range(self._mdsystem.CG.ncell):
            cell_ptr = &(self._mdsystem.CG.data[i])

            for j in range(cell_ptr.nr[0]):
                k = cell_ptr.data[j].id
                gradients[k*3  ] = cell_ptr.data[j].frc[0]
                gradients[k*3+1] = cell_ptr.data[j].frc[1]
                gradients[k*3+2] = cell_ptr.data[j].frc[2]

        ene = self._mdsystem.sysUV
        return ene, np.asarray(gradients)


def test_load(tname, cname):
    import parmed as pmd
    import numpy as np

    parm = pmd.load_file(tname, cname)
    natoms = len(parm.atoms) 

    #initialize the charge array
    charge_array = np.zeros(natoms)

    for i in range(natoms):
        charge_array[i]= parm.atoms[i].charge 

    #create the double* as a Memory view type
    cdef double[:] charge_view = charge_array

    cdef double[:] charge_results = np.zeros(natoms)

    CpyDVec(&charge_view[0],natoms,&charge_results[0])

    #for i in range(natoms):
    #    print(charge_results[i])

    #if we want we can return the memory slice as a dictionary
    #data = dict()
    #data["charges"] = np.asarray(charge_results)

    #otherwise we return it  as a numpy array
    data = np.asarray(charge_results)
    return data


if __name__ == '__main__':
    test_load('vAla3.prmtop', 'vAla3.rst7')
