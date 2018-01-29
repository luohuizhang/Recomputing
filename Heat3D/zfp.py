from __future__ import division
import sys
import array
import adios as ad
import os



if len(sys.argv) < 2 :
    print 'usage: python %s <file.bp>' % sys.argv[0]
    sys.exit(0)

statinfo=os.stat(sys.argv[1])
fsize=statinfo.st_size
print fsize
fin = open(sys.argv[1], 'rb')
vals = array.array('d',(int(fsize/16))*[0])

fin.readinto(vals)
fin.close()







print "\n>>> Writing ...\n"

ad.init_noxml()

ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 100);

fname_writer = sys.argv[1]+".bp"

fw = ad.writer(fname_writer)

fw.declare_group('group', method='POSIX1')


fw['dpot'] = vals
fin.close()


fw.close()

