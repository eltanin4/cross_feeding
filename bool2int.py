def bool2int(x):
    y = 0
    for i,j in enumerate(x):
        y += int(j) << i
    return y

def int2bool(x):
    y = [ int(e) for e in bin( x )[ 2: ][ :: -1 ] ]
    y += [0] * ( len( rxnMat ) - len( y ) )
    return np.array( y )
