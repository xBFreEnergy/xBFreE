"""

Author: Mario S. Valdés Tresanco
Copyright: 2023
License: MIT

"""
# import ctypes
#
# import numpy as np
# from mpi4py import MPI
#
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# name = MPI.Get_processor_name()
#
# # create new communicators for each node according to the rank type shared
# nodes_comm = comm.Split_type(MPI.COMM_TYPE_SHARED, key=rank)
# nodes_rank = nodes_comm.Get_rank()
#
# print('rank =', rank,  'nodes_comm.rank =',nodes_comm.rank, 'nodes_comm.Get_rank() =', nodes_comm.Get_rank(),
#       'name =', name)


# import numpy as np
# from mpi4py import MPI
# import random
#
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# root = 0
#
# local_array = [rank] * random.randint(2, 5)
# print("rank: {}, local_array: {}".format(rank, local_array))
#
# sendbuf = np.array(local_array)
#
# # Collect local array sizes using the high-level mpi4py gather
# sendcounts = np.array(comm.gather(len(sendbuf), root))
#
# if rank == root:
#     print("sendcounts: {}, total: {}".format(sendcounts, sum(sendcounts)))
#     recvbuf = np.empty(sum(sendcounts), dtype=int)
# else:
#     recvbuf = None
#
# comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=root)
# if rank == root:
#     print("Gathered array: {}".format(recvbuf))

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

nodes_comm = comm.Split_type(MPI.COMM_TYPE_SHARED, key=rank)
nodes_rank = nodes_comm.Get_rank()

if nodes_rank == 0:
    data = {'a': 7, 'b': 3.14}
else:
    data = None
data = comm.gather(data, root=0)
comm.bcast(data)

# comm.Barrier()

if rank == 0:
    # for i in range(size):
    #     assert data[i] == (i+1)**2
    print(f"{data = }, {rank = }")

print('####')

print(data)
