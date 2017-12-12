#===========================================================
#
#  PROJECT: vasp_unfold
#  FILE:    utils.py
#  AUTHOR:  Milan Tomic
#  EMAIL:   tomic@th.physik.uni-frankfurt.de
#  VERSION: 1.4
#  DATE:    December 12th 2017
#
#===========================================================

import numpy as np
import argparse
import sys
import fractions
import traceback

        
def post_error(error_info, show_traceback=False):
    '''Write error message to sderr and exit program
    '''
    message = '\nError: '+error_info+'\n\n'
    
    if show_traceback:
        # Show exception traceback
        tb_msg = traceback.format_exc()
        
        max_line_len = max(len(line) for line in tb_msg.split('\n'))
        
        message += 'Operation failed due to the following exception:\n'
        message += '='*max_line_len
        message += '\n' + traceback.format_exc()
        message += '='*max_line_len
    
    sys.stderr.write(message)
    
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
    and array of three fractions is returned.
    '''
    try:
        return [fractions.Fraction(s) if s.strip() != "0" else 1 for s in tstring.split(',')]
    except:
        post_error('Unable to parse string: "{0}". Check help for valid '
                   'translation generator specification'.format(tstring))

def version(vstring):
    '''Parse string containing version information.
    Valid form has dot separated digits. Returns tuple
    of integers which are to be lexicographically compared.
    '''
    try:
        return tuple(int(d) for d in vstring.split('.'))
    except:
        post_error('Unable to parse string: "{0}". The valid version '
                   'is composed of dot separated digists'.format(vstring))
                   
def lcm(a, b):
    '''Return lowest common multiple.'''
    return a * b // fractions.gcd(a, b)


def lcmm(*args):
    '''Return lcm of args.'''   
    return reduce(lcm, args)
    

def fraction_order(fract):
    '''Return the lowest integer multiple of a fraction 
    that yields an integer.
    '''
    if not isinstance(fract, fractions.Fraction):
        fract = fractions.Fraction(fract)
        
    return lcm(fract.numerator, fract.denominator) / fract.numerator
    

def frac_translation_order(trans):
    '''Returns the lowest integer multiple of a fractional 
    translation that yields translation with integer components.
    '''
    order = [fraction_order(ti) for ti in trans]
    
    return lcmm(*order)
    