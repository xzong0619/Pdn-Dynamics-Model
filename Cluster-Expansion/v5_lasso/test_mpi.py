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

'''
Commend:
mpirun -n 4 python test_MPI.py
'''

'''
Caviness output


1 of 4
rank 1  has data: 3.
Core 1 Before gather.
Core 1 After gather.
0 of 4
we will be scattering:[1, 2, 9, 64]
rank 0  has data: 2.
Core 0 Before gather.
Core 0 After gather.
master:[2, 3, 10, 65]
2 of 4
rank 2  has data: 10.
Core 2 Before gather.
Core 2 After gather.
3 of 4
rank 3  has data: 65.
Core 3 Before gather.
Core 3 After gather.


'''
'''
Mac output

0 of 4
we will be scattering:[1, 2, 9, 64]
1 of 4
rank 0  has data: 2.
Core 0 Before gather.
2 of 4
rank 2  has data: 10.
Core 2 Before gather.
Core 2 After gather.
rank 1  has data: 3.
Core 1 Before gather.
Core 1 After gather.
3 of 4
rank 3  has data: 65.
Core 3 Before gather.
Core 3 After gather.
Core 0 After gather.
master:[2, 3, 10, 65]
'''