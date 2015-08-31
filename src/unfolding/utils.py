#===========================================================
#
#  PROJECT: vasp_unfold
#  FILE:    utils.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.3
#  DATE:    Aug 31st 2015
#
#===========================================================

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
    
    def __init__(self, fname, comment=None):
        '''Constructor opens the file in the read mode'''
        super(Getlines, self).__init__(fname, 'r')
        self.comment = comment
    
    
    def __advance__(self, advance_func, eof_error):
        '''Advances along file line by line, skipping empty
        lines and stripping trailing and leading whitspaces 
        and comments. advance_func specifies function which
        provides contents of next line. If eof_error is True
        error will be posted and execution terminated if end
        of file is reached.
        '''
        line = ''
        
        while len(line) == 0:
            line = advance_func()
            
            if line == '':
                if eof_error is True:
                    post_error('Reached end of: {0}.'.format(self.name))
                else:
                    return None
            if self.comment is not None:
                # Strip comments
                line = line[:line.find(self.comment)]
            
            # Strip spaces
            line = line.strip()
            
        return line
        
        
    def next(self):
        '''Overides file.next() so that empty lines are skipped
        and leading and trailing whitespaces and comments stripped.
        '''
        return self.__advance__(super(Getlines, self).next, False)
    
        
    def readline(self, eof_error=True):
        '''Overides file.readline() so that empty lines are skipped
        and leading and trailing whitespaces and comments stripped. 
        Optionally, if eof_error is True, error will be posted and 
        execution terminated in case end of file is reached
        '''
        return self.__advance__(super(Getlines, self).readline, eof_error)

    
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
    
