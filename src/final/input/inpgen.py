#!/usr/bin/python
import sys
import os
import random

# global defines
rstart = -100
rend = 100

def main():
    if(len(sys.argv) != 2):
        print 'Usage: %s dim' %(sys.argv[0])
    else:
        random.seed(42)
        dim = int(sys.argv[1])
        inp_fname = 'test_input_%d.in' %(dim)
        inp = open(inp_fname, 'w')
        inp.write('%d\n' %(dim))
        for row in range(0, dim):
            for col in range(0, dim):
                inp.write('%d\t' %(random.randint(rstart, rend)))
            inp.write('\n')
        inp.close()
        print 'File: %s has been created!' %(inp_fname)

if __name__ == '__main__':
    main()
