from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
print(size)
rank = comm.Get_rank()
print(rank)
if rank == 0:
   data = [(x+1)**x for x in range(size)]
   print('we will be scattering:{}'.format(data))
else:
   data = None
   
data = comm.scatter(data, root=0)
data += 1
print('rank {}  has data: {}.'.format(rank, data))

newData = comm.gather(data,root=0)

if rank == 0:
   print('master:{}'.format(newData))