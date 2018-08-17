from __future__ import print_function
import os
import sys
import numpy as np

# desdb is not needed in all scenarios
try:
    import desdb
except ImportError:
    pass


def read_blacklist(fname):
    import esutil as eu
    dt=[('expnum','i8'),
        ('ccd','i8')]

    #print("reading:",fname)
    with eu.recfile.Recfile(fname, 'r', delim=' ', dtype=dt) as fobj:
        data=fobj.read()
    """
    dt=[('expnum','i8'),
        ('ccd','i8')]    
    data = np.genfromtxt(fname,dtype=dt)
    """
    return data

def read_blacklist_as_dict(fname):
    data=read_blacklist(fname)

    bigind=make_bigind(data['expnum'], data['ccd'])
    d={}
    for i in xrange(data.size):
        d[bigind[i]] = data[i]

    return d

def get_corrupted_blacklist():
    subdirs=['EXTRA','blacklists']
    dir=desdb.files.get_dir_generic(subdirs)

    fname=os.path.join(dir, 'corrupted-y1.txt')
    return read_blacklist_as_dict(fname)

def get_exp_blacklists():
    subdirs=['EXTRA','blacklists']
    dir=desdb.files.get_dir_generic(subdirs)

    # Eli's flags go to 2**9

    #print("reading blackists")
    ldict={}
    fname=os.path.join(dir, 'ghost-scatter-sv-uniq.txt')
    ldict['ghost-sv'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**10}

    fname=os.path.join(dir, 'ghost-scatter-y1-uniq.txt')
    ldict['ghost-y1'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**11}

    fname=os.path.join(dir, 'noise-y1-uniq.txt')
    ldict['noise-y1'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**12}

    fname=os.path.join(dir, 'streak-sv-uniq.txt')
    ldict['streak-sv'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**13}

    fname=os.path.join(dir, 'streak-y1-uniq.txt')
    ldict['streak-y1'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**14}

    return ldict

def remove_corrupted(srclist):
    blacklist = get_corrupted_blacklist()

    new_srclist=[]
    for s in srclist:
        if s['bigind'] not in blacklist:
            new_srclist.append( s )
        else:
            print("    found in corrupted blacklist")

    print("kept %d/%d after "
          "removing corrupted" % (len(new_srclist),len(srclist)))
    return new_srclist

def add_blacklist_flags(srclist):
    """
    bigind and flags must be present already
    """
    blacklists = get_exp_blacklists()
    for s in srclist:
        for bname,blist in blacklists.iteritems():
            if s['bigind'] in blist['blacklist']:
                print("    found in blacklist:",bname)
                s['flags'] |= blist['flag']

def add_bigind(srclist):
    for s in srclist:
        expname=s['expname']
        expnum=int( expname.split('_')[1] )
        s['bigind'] = make_bigind(expnum, s['ccd'])

def make_bigind(expnum, ccdnum):
    return expnum + ccdnum*10**7
