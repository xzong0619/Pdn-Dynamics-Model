from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print("%d of %d" % (comm.Get_rank(), comm.Get_size()))

if rank == 0:
   data = [(x+1)**x for x in range(size)]
   print('we will be scattering:{}'.format(data))
else:
   data = None
   
data = comm.scatter(data, root=0)
data += 1
print('rank {}  has data: {}.'.format(rank, data))

print('Core {} Before gather.'.format(rank))

newData = comm.gather(data,root=0)

print('Core {} After gather.'.format(rank))
if rank == 0:
   print('master:{}'.format(newData))