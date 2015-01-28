import numpy as np
import argparse
import sys

        
def post_error(error_info):
    '''Write error message to sderr and exit program
    '''
    sys.stderr.write('\nError: '+error_info+'\n\n')
    
    exit()
    
    
class Getlines(file):
    '''Small wrapper of the Python's built-in file
    class. It's purpose is to skip empty lines
    while reading the file'''
    
    def __init__(self, fname):
        '''Constructor opens the file in the read mode'''
        super(Getlines, self).__init__(fname, 'r')
    
    def get_next(self, eof_error=True):
        '''Advance to the next non-empty line and return
        it, while stripping leading and trailing white-
        spaces and the newline. When end of file is reached
        error will be posted and exectuion terminated if
        eof_error is set to True. Otherwise, None is returned.'''
        
        line = ''
        
        while len(line) == 0:
            line = self.readline()
            
            # If we reach end of file, we can
            # either post eror and abort execution
            # or just return None
            if line == '':
                if eof_error:
                    post_error('File "{0}" incomplete'.format(self.name))
                else:
                    return None
            else:
                line = line.strip()
                
        return line

    
def translation(tstring):
    '''Parse string describing fractional translations.
    Valid form has comma separated components, where 
    every component either 0 or 1/n, with n positive
    integer. During parsing 0 is replaced with 1/1
    and integer array containing denominators is returned.
    '''
    try:        
        tstring = tstring.replace('0', ' 1/1')

        return np.array([int(s.split('/')[1]) for s in tstring.split(',')])
    except:
        post_error('Unable to parse string: "{0}". Check help for valid '
                   'translation generator specification'.format(tstring))
    
    
def gcd(a, b):
    '''Return greatest common divisor using Euclid's Algorithm.'''
    while b:      
        a, b = b, a % b
    return a


def lcm(a, b):
    '''Return lowest common multiple.'''
    return a * b // gcd(a, b)


def lcmm(*args):
    '''Return lcm of args.'''   
    return reduce(lcm, args)
    
