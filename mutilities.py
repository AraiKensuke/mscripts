import multiprocessing as _mp


def nosharedDatMP(target, thrItr, maxusecpus=1):
    #  
    #  for running simple functions using multiple threads
    #    threads are differentiated by different values of items in thrItr
    #    that piece of information should be all that is necessary by the
    #    thread to function
    #  TODO:  how does the target interact with global variables
    #  
    if maxusecpus > 1:
        ncpus   = _mp.cpu_count()
        usecpus = ncpus if ncpus <= maxusecpus else maxusecpus
        print usecpus

        pool    = _mp.Pool(usecpus)    #  max of 10 cores
        pool.map(target, thrItr)
#        pool = []
    #     for itr in thrItr:   #  key is exptDate
    #         pool.append(_mp.Process(target=target, args=(itr,)))
    #     for p in pool:
    #         p.start()
    #     for p in pool:
    #         p.join()
    else:
        for itr in thrItr:   #  key is exptDate
            target(itr)

#def sharedDatMP(target, thrItr, sharedMemItems, maxusecpus=1):
#     if maxusecpus == 1:
#         ncpus   = _mp.cpu_count()
#         usecpus = ncpus if ncpus <= maxusecpus else maxusecpus
#         print "usecpus:   %d" % usecpus
#         pool = []

#         shalldat = _mp.RawArray(ctypes.c_short, size)
#         _alldat    = shmem_as_ndarray(shalldat).reshape(rows,cols)
#         _alldat[:, :] = alldat[:, :]

#         for exptDate in n2u:   #  key is exptDate
#             pool.append(_mp.Process(target=_filterAlldat, args=(exptDate,)))
#         for p in pool:
#             p.start()
#         for p in pool:
#             p.join()

#     else:
#         for exptDate in n2u:   #  key is exptDate
#             chi2ForExptDate(exptDate)
    
