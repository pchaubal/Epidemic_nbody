import time

def timeit(method):
    def timed(*args, **kw):
        print ('Running %r' % method.__name__)
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print ('%r %2.2f sec' % \
              (method.__name__, te-ts))
        return result
    return timed
