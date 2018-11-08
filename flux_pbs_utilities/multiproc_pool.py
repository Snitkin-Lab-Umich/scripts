import multiprocessing as mp

def list_append(count, id, out_list):
        """
        Creates an empty list and then appends a
        random number to the list 'count' number
        of times. A CPU-heavy operation!
        """
	for i in range(count):
                out_list.append(random.random())

def run_command(i):
    print "Running: %s" % i
    #os.system(i)
    done = "done: %s" % i
    return done

num_workers = mp.cpu_count()  
size = 10000000   # Number of random numbers to add
#procs = 12
pool = mp.Pool(num_workers)
for i in range(0, num_workers):
    out_list = list()    
    pool.apply_async(target=run_command, args=(i,))

pool.close()
pool.join()
