from __future__ import print_function
import os 
import pymp
print(os.cpu_count()) 
ex_array = pymp.shared.array(100)
for i in range(100):
    ex_array[i]=i
with pymp.Parallel(100) as p:
    for index in p.range(0, 100):
        ex_array[index] = ex_array[index]*10
        # The parallel print function takes care of asynchronous output.
        #p.print('Yay! {} done!'.format(index))
        p.print(p.thread_num)
print(ex_array)        