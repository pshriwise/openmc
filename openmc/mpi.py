try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except ImportError:
    from unittest.mock import Mock
    MPI = Mock()
    from openmc.dummy_comm import DummyCommunicator
    comm = DummyCommunicator()


def with_barrier(rank=0):
    """
    Decorator indicating that a function should always run
    on a specified rank with a barrier for all other procs

    func : Callable
        The function to be called
    rank : int
        Process rank that will do the work (default 0)
    """
    def dec(func):
        def wrapper(*args, **kwargs):
            if comm.rank == rank:
                func(*args, **kwargs)
            comm.barrier()
        return wrapper
    return dec