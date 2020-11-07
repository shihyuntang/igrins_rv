import numpy as np

def sigmaclip(l_return,l_ref,n):
    go = True; l = l_return; lm = l_ref;
    while go == True:
        lenbefore = len(lm)
        m = np.mean(lm)
        s = np.std(lm)
        l = l[(lm < m+n*s) & (lm > m-n*s)]
        lm = lm[(lm < m+n*s) & (lm > m-n*s)]
        if lenbefore == len(lm):
            return l

def sigmaclip_above(l_return,l_ref,n):
    go = True; l_out = np.array([]); l_base = np.array([]);
    while go == True:
        lenbefore = len(l_ref);
        for clic in range(15):
            l_refchunk = l_ref[np.int(clic*len(l_ref)/15.):np.int((clic+1)*len(l_ref)/15.)]
            l_returnchunk = l_return[np.int(clic*len(l_ref)/15.):np.int((clic+1)*len(l_ref)/15.)]
            m = np.median(l_refchunk);     s = np.std(l_refchunk);
            l_returnchunk = l_returnchunk[(l_refchunk < m + n*s)]; l_refchunk = l_refchunk[(l_refchunk < m + n*s)];
            l_base = np.concatenate((l_base,l_refchunk));  l_out = np.concatenate((l_out,l_returnchunk));
        if lenbefore == len(l_base):
            return l_out
        l_return = l_out; l_ref = l_base;
        l_out = np.array([]); l_base = np.array([]);

def basicclip_above(l_return,l_ref,nslices):
    n = 3
    s2nd  = np.sort(l_ref)
    s2ndd = s2nd[::-1]
    tops  = s2ndd[n-1] # third flux
    l_return =  l_return[(l_ref < tops)]
    l_ref    =  l_ref[   (l_ref < tops)]

    go = True; l_out = np.array([]); l_base = np.array([]);
    for clic in range(nslices):
        l_refchunk    = l_ref[   np.int(clic*len(l_ref)/np.float(nslices)) : np.int((clic+1)*len(l_ref)/np.float(nslices))]
        l_returnchunk = l_return[np.int(clic*len(l_ref)/np.float(nslices)) : np.int((clic+1)*len(l_ref)/np.float(nslices))]
        #xj = np.arange(len(l_refchunk))
        #q = np.polyfit(xj,l_refchunk,2)
        #f = np.poly1d(q)
        #l_returnchunk = l_returnchunk[(l_refchunk/f(xj) < max(l_refchunk/f(xj)))];
        l_returnchunk = l_returnchunk[(l_refchunk < np.max(l_refchunk))];
        l_out = np.concatenate((l_out,l_returnchunk));
    return l_out


def singler(xin):
    # Index selector for use during isolation of CRs and misfits after first run at optimization
    xout = np.array([])
    for x in xin:
        if (x+1 in xin and x-1 in xin) or (x+1 in xin and x+2 in xin) or  (x-1 in xin and x-2 in xin):
            continue
        else:
            xout = np.append(xout, np.int(x))
    return np.array(xout,dtype=int)

def multipler(xin):
    # Index selector for use during isolation of CRs and misfits after first run at optimization
    xout = np.array([])
    for x in xin:
        if (x+1 in xin and x-1 in xin) or (x+1 in xin and x+2 in xin) or  (x-1 in xin and x-2 in xin):
            xout = np.append(xout, np.int(x))
    return np.array(xout,dtype=int)
