#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Number theoretic functions
"""

from __future__ import unicode_literals
from __future__ import absolute_import
import six
from six.moves import map
from six.moves import range

import sympy
from itertools import combinations

from mathics.builtin.base import Builtin, Test
from mathics.core.expression import (
    Expression, Integer, Rational, Symbol, from_python)
from mathics.core.numbers import (convert_base, convert_int_to_digit_list)
import mpmath
import string
import re
from mathics.builtin.algebra import Denominator
from django.utils.regex_helper import contains
class PowerMod(Builtin):
    """
    <dl>
    <dt>'PowerMod[$x$, $y$, $m$]'
        <dd>computes $x$^$y$ modulo $m$.
    </dl>

    >> PowerMod[2, 10000000, 3]
     = 1
    >> PowerMod[3, -2, 10]
     = 9
    >> PowerMod[0, -1, 2]
     : 0 is not invertible modulo 2.
     = PowerMod[0, -1, 2]
    >> PowerMod[5, 2, 0]
     : The argument 0 should be nonzero.
     = PowerMod[5, 2, 0]

    'PowerMod' does not support rational coefficients (roots) yet.
    """

    attributes = ('Listable',)

    messages = {
        'ninv': "`1` is not invertible modulo `2`.",
    }

    def apply(self, a, b, m, evaluation):
        'PowerMod[a_Integer, b_Integer, m_Integer]'

        a_int = a
        m_int = m
        a, b, m = a.value, b.value, m.value
        if m == 0:
            evaluation.message('PowerMod', 'divz', m)
            return
        if b < 0:
            b = -b
            try:
                a = int(sympy.invert(a, m))
            except sympy.polys.polyerrors.NotInvertible:
                evaluation.message('PowerMod', 'ninv', a_int, m_int)
                return
        return Integer(pow(a, b, m))


class Mod(Builtin):
    """
    <dl>
    <dt>'Mod[$x$, $m$]'
        <dd>returns $x$ modulo $m$.
    </dl>

    >> Mod[14, 6]
     = 2
    >> Mod[-3, 4]
     = 1
    >> Mod[-3, -4]
     = -3
    >> Mod[5, 0]
     : The argument 0 should be nonzero.
     = Mod[5, 0]
    """

    attributes = ('Listable', 'NumericFunction')

    def apply(self, n, m, evaluation):
        'Mod[n_Integer, m_Integer]'

        n, m = n.get_int_value(), m.get_int_value()
        if m == 0:
            evaluation.message('Mod', 'divz', m)
            return
        return Integer(n % m)


class EvenQ(Test):
    """
    <dl>
    <dt>'EvenQ[$x$]'
        <dd>returns 'True' if $x$ is even, and 'False' otherwise.
    </dl>

    >> EvenQ[4]
     = True
    >> EvenQ[-3]
     = False
    >> EvenQ[n]
     = False
    """

    def test(self, n):
        value = n.get_int_value()
        return value is not None and value % 2 == 0


class OddQ(Test):
    """
    <dl>
    <dt>'OddQ[$x$]'
        <dd>returns 'True' if $x$ is odd, and 'False' otherwise.
    </dl>

    >> OddQ[-3]
     = True
    >> OddQ[0]
     = False
    """

    def test(self, n):
        value = n.get_int_value()
        return value is not None and value % 2 != 0


class GCD(Builtin):
    """
    <dl>
    <dt>'GCD[$n1$, $n2$, ...]'
        <dd>computes the greatest common divisor of the given integers.
    </dl>

    >> GCD[20, 30]
     = 10
    >> GCD[10, y]
     = GCD[10, y]

    'GCD' is 'Listable':
    >> GCD[4, {10, 11, 12, 13, 14}]
     = {2, 1, 4, 1, 2}

    'GCD' does not work for rational numbers and Gaussian integers yet.
    """

    attributes = ('Flat', 'Listable', 'OneIdentity', 'Orderless')

    def apply(self, ns, evaluation):
        'GCD[ns___Integer]'

        ns = ns.get_sequence()
        result = 0
        for n in ns:
            value = n.get_int_value()
            if value is None:
                return
            result = sympy.gcd(result, value)
        return Integer(result)

# FIXME: Previosuly this used gmpy's gcdext. sympy's gcdex is not as powerful
# class ExtendedGCD(Builtin):
#    """
#    >> ExtendedGCD[10, 15]
#     = {5, {-1, 1}}
#
#    'ExtendedGCD' works with any number of arguments:
#    >> ExtendedGCD[10, 15, 7]
#     = {1, {-3, 3, -2}}
#
#    Compute the greated common divisor and check the result:
#    >> numbers = {10, 20, 14};
#    >> {gcd, factors} = ExtendedGCD[Sequence @@ numbers]
#     = {2, {3, 0, -2}}
#    >> Plus @@ (numbers * factors)
#     = 2
#
#    'ExtendedGCD' does not work for rational numbers and Gaussian integers yet
#    """
#
#    attributes = ('Listable',)
#
#    def apply(self, ns, evaluation):
#        'ExtendedGCD[ns___Integer]'
#
#        ns = ns.get_sequence()
#        result = 0
#        coeff = []
#        for n in ns:
#            value = n.get_int_value()
#            if value is None:
#                return
#            new_result, c1, c2 = sympy.gcdex(result, value)
#            result = new_result
#            coeff = [c * c1 for c in coeff] + [c2]
#            return Expression('List', Integer(result), Expression(
#                'List', *(Integer(c) for c in coeff)))


class LCM(Builtin):
    """
    <dl>
    <dt>'LCM[$n1$, $n2$, ...]'
        <dd>computes the least common multiple of the given integers.
    </dl>

    >> LCM[15, 20]
     = 60
    >> LCM[20, 30, 40, 50]
     = 600
    """

    attributes = ('Flat', 'Listable', 'OneIdentity', 'Orderless')

    def apply(self, ns, evaluation):
        'LCM[ns___Integer]'

        ns = ns.get_sequence()
        result = 1
        for n in ns:
            value = n.get_int_value()
            if value is None:
                return
            result = sympy.lcm(result, value)
        return Integer(result)


class FactorInteger(Builtin):
    """
    <dl>
    <dt>'FactorInteger[$n$]'
        <dd>returns the factorization of $n$ as a list of factors and exponents.
    </dl>

    >> factors = FactorInteger[2010]
     = {{2, 1}, {3, 1}, {5, 1}, {67, 1}}
    To get back the original number:
    >> Times @@ Power @@@ factors
     = 2010
    'FactorInteger' factors rationals using negative exponents:
    >> FactorInteger[2010 / 2011]
     = {{2, 1}, {3, 1}, {5, 1}, {67, 1}, {2011, -1}}
    """

    # TODO: GausianIntegers option
    # e.g. FactorInteger[5, GaussianIntegers -> True]

    def apply(self, n, evaluation):
        'FactorInteger[n_]'

        if isinstance(n, Integer):
            factors = sympy.factorint(n.value)
            factors = sorted(six.iteritems(factors))
            return Expression('List', *(Expression('List', factor, exp)
                                        for factor, exp in factors))

        elif isinstance(n, Rational):
            factors, factors_denom = list(map(
                sympy.factorint, n.value.as_numer_denom()))
            for factor, exp in six.iteritems(factors_denom):
                factors[factor] = factors.get(factor, 0) - exp
            factors = sorted(six.iteritems(factors))
            return Expression('List', *(Expression('List', factor, exp)
                                        for factor, exp in factors))
        else:
            return evaluation.message('FactorInteger', 'exact', n)


class IntegerExponent(Builtin):
    """
    <dl>
    <dt>'IntegerExponent[$n$, $b$]'
        <dd>gives the highest exponent of $b$ that divides $n$.
    </dl>

    >> IntegerExponent[16, 2]
     = 4

    >> IntegerExponent[-510000]
     = 4

    >> IntegerExponent[10, b]
     = IntegerExponent[10, b]
    """

    rules = {
        'IntegerExponent[n_]': 'IntegerExponent[n, 10]',
    }

    messages = {
        'int': 'Integer expected at position 1 in `1`',
        'ibase': 'Base `1` is not an integer greater than 1.',
    }

    def apply(self, n, b, evaluation):
        'IntegerExponent[n_Integer, b_Integer]'

        py_n, py_b = n.to_python(), b.to_python()
        expr = Expression('IntegerExponent', n, b)

        if not isinstance(py_n, six.integer_types):
            evaluation.message('IntegerExponent', 'int', expr)
        py_n = abs(py_n)

        if not (isinstance(py_b, int) and py_b > 1):
            evaluation.message('IntegerExponent', 'ibase', b)

        # TODO: Optimise this (dont need to calc. base^result)
        # NOTE: IntegerExponent[a,b] causes a Python error here when a or b are
        # symbols
        result = 1
        while py_n % (py_b ** result) == 0:
            result += 1

        return from_python(result - 1)

def check_finite_decimal(denominator):
    # The rational number is finite decimal if the denominator has form 2^a * 5^b
    while denominator % 5 == 0:
        denominator = denominator / 5
    
    while denominator % 2 == 0:
        denominator = denominator / 2
            
    return True if denominator == 1 else False

def convert_repeating_decimal(numerator, denominator, base):
    head = [x for x in str(numerator // denominator)]
    tails = []
    subresults = [numerator % denominator]
    numerator %= denominator
    
    while numerator != 0:  # only rational input can go to this case
        numerator *= base
        result_digit, numerator = divmod(numerator, denominator)
        tails.append(str(result_digit))
        if numerator not in subresults:
            subresults.append(numerator)           
        else:
            break
        
    for i in range(len(head) - 1, -1, -1):
        j = len(tails) - 1
        if head[i] != tails[j]:
            break;
        else:
            del(tails[j])
            tails.insert(0, head[i])
            del(head[i])
            j = j - 1
    
    # truncate all leading 0's        
    if all(elem == '0' for elem in head):
        for i in range(0, len(tails)):
            if tails[0] == '0':
               tails = tails[1:] + [str(0)]
            else:
               break 
    return (head, tails)

def convert_float_base(x, base, precision=10):
    
    length_of_int = 0 if x == 0 else int(mpmath.log(x, base))
    iexps = list(range(length_of_int, -1, -1))
    
    def convert_int(x, base, exponents):
        out = []
        for e in range(0, exponents + 1):
            d = x % base
            out.append(d)
            x = x / base
            if x == 0:
                break
        out.reverse()
        return out

    def convert_float(x, base, exponents):
        out = []
        for e in range(0, exponents):
            d = int(x * base)
            out.append(d)
            x = (x * base) - d
            if x == 0:
                break
        return out
        # print(out,out.reverse())
    int_part = convert_int(int(x), base, length_of_int)
    if isinstance(x, (float, sympy.Float)):
        fexps = list(range(-1, -int(precision + 1), -1))
        # real_part = convert(x - int(x), base, fexps, True, int_part)
        real_part = convert_float(x - int(x), base, precision + 1)
        return int_part + real_part
    elif isinstance(x, six.integer_types):
        return int_part
    else:
        raise TypeError(x)

class RealDigits(Builtin):
    """
    <dl>
    <dt>'RealDigits[$n$]'
        <dd>finds a list of the digits in the approximate real number $n$, together with the number of digits that are to the left of the decimal point.
    
    <dt>'RealDigits[$n$, $b$]'
        <dd>finds a list of base‐$b$ digits in $n$.
    
    <dt>'RealDigits[$n$, $b$, $len$]'
        <dd>finds a list of $len$ digits.
    
    <dt>'RealDigits[$n$, $b$, $len$, $p$]'
        <dd>finds $len$ digits starting with the coefficient of $b$$p$
    </dl>

    >> RealDigits[123.55555]
     = {{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}
    
    >> RealDigits[0.000012355555]
     = {{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, -4}
    
    >> RealDigits[-123.55555]
     = {{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}
    
    #> RealDigits[0.004]
     = {{4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, -2}
     
    #> RealDigits[-1.25, -1]
     : Base -1 is not a real number greater than 1.
     = RealDigits[-1.25, -1] 
     
    #> RealDigits[Pi, 10, 25]
     = {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, 1}
    
    #> RealDigits[19 / 7, 10, 25]
     = {{2, 7, 1, 4, 2, 8, 5, 7, 1, 4, 2, 8, 5, 7, 1, 4, 2, 8, 5, 7, 1, 4, 2, 8, 5}, 1} 
     
    #> RealDigits[19 / 7]
     = {{2, {7, 1, 4, 2, 8, 5}}, 1}
     
    #> RealDigits[100 / 21]
     = {{{4, 7, 6, 1, 9, 0}}, 1}
     
    #> RealDigits[1.234, 2, 15]
     = {{1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1}, 1}
     
    #> RealDigits[Pi, 10, 20, -5]
     = {{9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, -4}
     
    #> RealDigits[Pi, 10, 20, 5]
     = {{0, 0, 0, 0, 0, 3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9}, 6}
     
    #> RealDigits[Pi, 10, 1, -10000]
     = {{8}, -9999}
     
    #> RealDigits[Pi]
     : The number of digits to return cannot be determined.
     = RealDigits[Pi]
    
    #> RealDigits[20 / 3]
     = {{{6}}, 1}
     
    #> RealDigits[3 / 4]
     = {{7, 5}, 0}
     
    #> RealDigits[23 / 4]
     = {{5, 7, 5}, 1}
     
    #> RealDigits[3 + 4 I]
     : The value 3 + 4 I is not a real number.
     = RealDigits[3 + 4 I]
     
    #> RealDigits[abc]
     = RealDigits[abc]
     
    #> RealDigits[abc, 2]
     = RealDigits[abc, 2]  
     
    #> RealDigits[45]
     = {{4, 5}, 2}
     
    #> RealDigits[{3.14, 4.5}]
     = {{{3, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1}, {{4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1}}
    
    #> RealDigits[123.45, 40]
     = {{3, 3, 18, 0, 0, 0, 0, 0, 0, 0}, 2}
     
    #> RealDigits[0.00012345, 2]
     = {{1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0}, -12}
    
    #> RealDigits[12345, 2, 4]
     = {{1, 1, 0, 0}, 14}
     
    #> RealDigits[123.45, 2, 15]
     = {{1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1}, 7}
     
    #> RealDigits[123.45, 10, 18]
     = {{1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Indeterminate, Indeterminate}, 3}
     
    #> RealDigits[0.000012345, 2]
     = {{1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1}, -16}
    
    #> RealDigits[3.14, 10, 1.5]
     : Non-negative machine-sized integer expected at position 3 in RealDigits[3.14, 10, 1.5].
     = RealDigits[3.14, 10, 1.5]
     
    #> RealDigits[3.14, 10, 1, 1.5]
     : Machine-sized integer expected at position 4 in RealDigits[3.14, 10, 1, 1.5].
     = RealDigits[3.14, 10, 1, 1.5]
     
    #> RealDigits[Pi, 10, 20, -5]
     = {{9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, -4}
     
    #> RealDigits[305.0123, 10, 17, 0]
     = {{5, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, Indeterminate, Indeterminate, Indeterminate}, 1}
     
    #> RealDigits[220, 140]
     = {{1, 80}, 2}
     
    #> RealDigits[Sqrt[3], 10, 50]
     = {{1, 7, 3, 2, 0, 5, 0, 8, 0, 7, 5, 6, 8, 8, 7, 7, 2, 9, 3, 5, 2, 7, 4, 4, 6, 3, 4, 1, 5, 0, 5, 8, 7, 2, 3, 6, 6, 9, 4, 2, 8, 0, 5, 2, 5, 3, 8, 1, 0, 3}, 1}
    
    #> RealDigits[0]
     = {{0}, 1}
    
    #> RealDigits[1]
     = {{1}, 1}
     
    #> RealDigits[0, 10, 5]
     = {{0, 0, 0, 0, 0}, 0}
     
    #> RealDigits[11/23]
     = {{{4, 7, 8, 2, 6, 0, 8, 6, 9, 5, 6, 5, 2, 1, 7, 3, 9, 1, 3, 0, 4, 3}}, 0}
    
    #> RealDigits[1/97]
     = {{{1, 0, 3, 0, 9, 2, 7, 8, 3, 5, 0, 5, 1, 5, 4, 6, 3, 9, 1, 7, 5, 2, 5, 7, 7, 3, 1, 9, 5, 8, 7, 6, 2, 8, 8, 6, 5, 9, 7, 9, 3, 8, 1, 4, 4, 3, 2, 9, 8, 9, 6, 9, 0, 7, 2, 1, 6, 4, 9, 4, 8, 4, 5, 3, 6, 0, 8, 2, 4, 7, 4, 2, 2, 6, 8, 0, 4, 1, 2, 3, 7, 1, 1, 3, 4, 0, 2, 0, 6, 1, 8, 5, 5, 6, 7, 0}}, -1}
    
    #> RealDigits[1/97, 2]
     = {{{1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}}, -6}
   
    #> RealDigits[1/197, 260, 5]
     = {{1, 83, 38, 71, 69}, 0}
     
    #> RealDigits[1/197, 260, 5, -6]
     = {{246, 208, 137, 67, 80}, -5}
     
    #> RealDigits[Pi, 260, 20]
     = {{3, 36, 211, 172, 124, 173, 210, 42, 162, 76, 23, 206, 122, 187, 23, 245, 241, 225, 254, 98}, 1}
  
    #> RealDigits[Pi, 260, 5]
     = {{3, 36, 211, 172, 124}, 1}
     
    #> RealDigits[1/3]
     = {{{3}}, 0}
     
    #> RealDigits[1/2, 7]
     = {{{3}}, 0}
     
    #> RealDigits[3/2, 7]
     = {{1, {3}}, 1}
    
    #> RealDigits[-3/2, 7]
     = {{1, {3}}, 1}
     
    #> RealDigits[3/2, 6]
     = {{1, 3}, 1}
     
    #> RealDigits[1, 7, 5]
     = {{1, 0, 0, 0, 0}, 1}
     
    #> RealDigits[I, 7]
     : The value I is not a real number.
     = RealDigits[I, 7]
     
    #> RealDigits[-Pi]
     : The number of digits to return cannot be determined.
     = RealDigits[-Pi]
     
    #> RealDigits[Round[x + y]]
     = RealDigits[Round[x + y]]
    
     
    """
    attributes = ('Listable',)
    
    messages = {
        'realx': 'The value `1` is not a real number.',
        'ndig': 'The number of digits to return cannot be determined.',
        'rbase': 'Base `1` is not a real number greater than 1.',
        'intnm': 'Non-negative machine-sized integer expected at position 3 in `1`.',
        'intm' : 'Machine-sized integer expected at position 4 in `1`.',
    }
    
    def apply_complex(self, n, var, evaluation):
        'RealDigits[n_Complex, var___]'
        return evaluation.message('RealDigits', 'realx', n)
    
    def apply_rational_with_base(self, n, b, evaluation):
        'RealDigits[n_Rational, b_Integer]'
        
        expr = Expression('RealDigits', n)
           
        py_n = abs(n.value)
        py_b = b.to_python()
        
        if check_finite_decimal(n.denominator().to_python()) == True and not py_b % 2 :
            return self.apply_2(n, b, evaluation)
        else:
            exp = int(mpmath.ceil(mpmath.log(py_n, py_b)))
            (head, tails) = convert_repeating_decimal(py_n.as_numer_denom()[0], py_n.as_numer_denom()[1], py_b)
            
            list_str = Expression('List')
            for x in head:
                if not x == '0':
                    list_str.leaves.append(from_python(int(x)))
            list_str.leaves.append(from_python(tails))
            
        return Expression('List', list_str, exp)
    
    def apply_rational_without_base(self, n, evaluation):
        'RealDigits[n_Rational]'
        
        return self.apply_rational_with_base(n, from_python(10), evaluation)
 
    def apply(self, n, evaluation):
        'RealDigits[n_]'
        
        # Handling the testcase that throw the error message and return the ouput that doesn't include `base` argument
        if isinstance(n, Symbol):
            return evaluation.message('RealDigits', 'ndig', n) if n.name.startswith('System`')  else Expression('RealDigits', n)
        
        if n.is_numeric():
            return self.apply_2(n, from_python(10), evaluation)
        
    def apply_2(self, n, b, evaluation, nr_elements=None, pos=None):
        'RealDigits[n_?NumericQ, b_Integer]'
        
        rational_no = True if isinstance(n, Rational) else False
        expr = Expression('RealDigits', n)
        py_b = b.to_python()
        if isinstance(n, (Expression, Symbol, Rational)):
            pos_len = abs(pos) + 1 if pos is not None and pos < 0 else 1
            if nr_elements is not None:
                n = Expression('N', n, int(mpmath.log(py_b ** (nr_elements + pos_len), 10)) + 1).evaluate(evaluation)
            else:
                if isinstance(n, Rational):
                    n = Expression('N', n).evaluate(evaluation)
                else:
                    return evaluation.message('RealDigits', 'ndig', expr)
        py_n = abs(n.value)
        
        if not py_b > 1:
            return evaluation.message('RealDigits', 'rbase', py_b)
                                      
        if isinstance(py_n, complex):
            return evaluation.message('RealDigits', 'realx', expr)
        
        if isinstance(n, Integer):
            display_len = int(mpmath.floor(mpmath.log(py_n, py_b))) if py_n != 0 and py_n != 1 else int(1)
        else:
            display_len = int(Expression('N', Expression('Round', Expression('Divide', Expression('Precision', py_n), Expression('Log', 10, py_b)))).evaluate(evaluation).to_python())
        
        exp = int(mpmath.ceil(mpmath.log(py_n, py_b))) if py_n != 0 and py_n != 1 else int(1)
        
        if py_n == 0 and nr_elements is not None:
            exp = 0
            
        digits = [] 
        
        if not py_b == 10: 
            digits = convert_float_base(py_n, py_b, display_len - exp)
            # truncate all the leading 0's
            for i in range(0, len(digits)):
                if(digits[0] != 0):
                    break
                else:
                    digits = digits[1:]
                    
            if not isinstance(n, Integer):
                if len(digits) > display_len:
                    digits = digits[:display_len-1]    
        else:
            # drop any leading zeroes
            for x in str(py_n):
                if x != '.' and (digits or x != '0'):
                    digits.append(x)
        if pos is not None:
            temp = exp
            exp = pos + 1
            move = temp - 1 - pos
            if move <= 0:
                digits = [0] * abs(move) + digits
            else:
                digits = digits[abs(move):]
                display_len = display_len - move 
                
        list_str = Expression('List')
        
        for x in digits:
            if x == 'e' or x == 'E':
                break 
            # Convert to Mathics' list format
            list_str.leaves.append(from_python(int(x)))  
        
        if rational_no is not True:
            while len(list_str.leaves) < display_len:
                list_str.leaves.append(from_python(0))
        
        if nr_elements is not None:
                    
            # display_len = nr_elements
            if len(list_str.leaves) >= nr_elements:
                # Truncate, preserving the digits on the right
                list_str = list_str.leaves[:nr_elements]   
            else:
                if isinstance(n, Integer):
                    while(len(list_str.leaves) < nr_elements):
                        list_str.leaves.append(from_python(0))
                else:
                    # Adding Indeterminate if the lenght is greater than the precision
                    while(len(list_str.leaves) < nr_elements):
                        list_str.leaves.append(from_python(Symbol('Indeterminate')))

        return Expression('List', list_str, exp)
    
    def apply_3(self, n, b, length, evaluation, pos=None):
        'RealDigits[n_?NumericQ, b_Integer, length_]'
        
        expr = Expression('RealDigits', n, b, length)
        
        if pos is not None: 
            expr.leaves.append(from_python(pos))
            
        if not(isinstance(length, Integer) and length.get_int_value() >= 0):
            return evaluation.message('RealDigits', 'intnm', expr)
        
        return self.apply_2(n, b, evaluation, nr_elements=length.get_int_value(), pos=pos)
    
    def apply_4(self, n, b, length, p, evaluation):
        'RealDigits[n_?NumericQ, b_Integer, length_, p_]'
        
        if not isinstance(p, Integer):
            return evaluation.message('RealDigits', 'intm', Expression('RealDigits', n, b, length, p))
        
        return self.apply_3(n, b, length, evaluation, pos=p.get_int_value())

class Prime(Builtin):
    """
    <dl>
    <dt>'Prime[$n$]'
        <dd>returns the $n$th prime number.
    </dl>

    >> Prime[1]
     = 2

    >> Prime[167]
     = 991
    """

    messages = {
        'intpp': 'Positive integer argument expected in `1`.',
    }

    def apply(self, n, evaluation):
        'Prime[n_]'
        n_int = n.to_python()
        if isinstance(n_int, int) and n_int > 0:
            return Integer(sympy.prime(n_int))

        expr = Expression('Prime', n)
        evaluation.message('Prime', 'intpp', expr)
        return


class PrimeQ(Builtin):
    """
    <dl>
    <dt>'PrimeQ[$n$]'
        <dd>returns 'True' if $n$ is a prime number.
    </dl>

    For very large numbers, 'PrimeQ' uses probabilistic prime testing, so it might be wrong sometimes
    (a number might be composite even though 'PrimeQ' says it is prime).
    The algorithm might be changed in the future.

    >> PrimeQ[2]
     = True
    >> PrimeQ[-3]
     = True
    >> PrimeQ[137]
     = True
    >> PrimeQ[2 ^ 127 - 1]
     = True

    #> PrimeQ[1]
     = False
    #> PrimeQ[2 ^ 255 - 1]
     = False

    All prime numbers between 1 and 100:
    >> Select[Range[100], PrimeQ]
     = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97}

    'PrimeQ' has attribute 'Listable':
    >> PrimeQ[Range[20]]
     = {False, True, True, False, True, False, True, False, False, False, True, False, True, False, False, False, True, False, True, False}
    """

    attributes = ('Listable',)

    def apply(self, n, evaluation):
        'PrimeQ[n_]'

        n = n.get_int_value()
        if n is None:
            return Symbol('False')

        n = abs(n)
        if sympy.isprime(n):
            return Symbol('True')
        else:
            return Symbol('False')


class CoprimeQ(Builtin):
    """
    <dl>
    <dt>'CoprimeQ[$x$, $y$]'
        <dd>tests whether $x$ and $y$ are coprime by computing their
        greatest common divisor.
    </dl>

    >> CoprimeQ[7, 9]
     = True

    >> CoprimeQ[-4, 9]
     = True

    >> CoprimeQ[12, 15]
     = False

    CoprimeQ also works for complex numbers
    >> CoprimeQ[1+2I, 1-I]
     = True

    >> CoprimeQ[4+2I, 6+3I]
     = False

    >> CoprimeQ[2, 3, 5]
     = True

    >> CoprimeQ[2, 4, 5]
     = False
    """
    attributes = ('Listable',)

    def apply(self, args, evaluation):
        'CoprimeQ[args__]'

        py_args = [arg.to_python() for arg in args.get_sequence()]
        if not all(isinstance(i, int) or isinstance(i, complex)
                   for i in py_args):
            return Symbol('False')

        if all(sympy.gcd(n, m) == 1 for (n, m) in combinations(py_args, 2)):
            return Symbol('True')
        else:
            return Symbol('False')


class PrimePowerQ(Builtin):
    """
    <dl>
    <dt>'PrimePowerQ[$n$]'
        <dd>returns 'True' if $n$ is a power of a prime number.
    </dl>

    >> PrimePowerQ[9]
     = True

    >> PrimePowerQ[52142]
     = False

    >> PrimePowerQ[-8]
     = True

    >> PrimePowerQ[371293]
     = True

    #> PrimePowerQ[1]
     = False
    """

    rules = {
        'PrimePowerQ[1]': 'False',
    }

    attributes = ('Listable', 'Protected', 'ReadProtected')

    # TODO: GaussianIntegers option
    """
    #> PrimePowerQ[5, GaussianIntegers -> True]
     = False
    """

    # TODO: Complex args
    """
    #> PrimePowerQ[{3 + I, 3 - 2 I, 3 + 4 I, 9 + 7 I}]
     = {False, True, True, False}
    """

    # TODO: Gaussian rationals
    """
    #> PrimePowerQ[2/125 - 11 I/125]
     = True
    """

    def apply(self, n, evaluation):
        'PrimePowerQ[n_]'
        n = n.get_int_value()
        if n is None:
            return Symbol('False')

        n = abs(n)
        if len(sympy.factorint(n)) == 1:
            return Symbol('True')
        else:
            return Symbol('False')


class PrimePi(Builtin):
    """
    <dl>
    <dt>'PrimePi[$x$]'
        <dd>gives the number of primes less than or equal to $x$.
    </dl>

    >> PrimePi[100]
     = 25

    >> PrimePi[-1]
     = 0

    >> PrimePi[3.5]
     = 2

    >> PrimePi[E]
     = 1
    """

    # TODO: Traditional Form

    def apply(self, n, evaluation):
        'PrimePi[n_?NumericQ]'
        result = sympy.ntheory.primepi(n.to_python(n_evaluation=evaluation))
        return from_python(result)


class NextPrime(Builtin):
    """
    <dl>
    <dt>'NextPrime[$n$]'
        <dd>gives the next prime after $n$.
    <dt>'NextPrime[$n$,$k$]'
        <dd>gives the $k$th  prime after $n$.
    </dl>

    >> NextPrime[10000]
     = 10007

    >> NextPrime[100, -5]
     = 73

    >> NextPrime[10, -5]
    = -2

    >> NextPrime[100, 5]
     = 113

    >> NextPrime[5.5, 100]
     = 563

    >> NextPrime[5, 10.5]
     = NextPrime[5, 10.5]
    """

    rules = {
        'NextPrime[n_]': 'NextPrime[n, 1]',
    }

    def apply(self, n, k, evaluation):
        'NextPrime[n_?NumericQ, k_?IntegerQ]'
        py_k = k.to_python(n_evaluation=evaluation)
        py_n = n.to_python(n_evaluation=evaluation)

        if py_k >= 0:
            return from_python(sympy.ntheory.nextprime(py_n, py_k))

        # Hack to get earlier primes
        result = n.to_python()
        for i in range(-py_k):
            try:
                result = sympy.ntheory.prevprime(result)
            except ValueError:
                # No earlier primes
                return from_python(-1 * sympy.ntheory.nextprime(0, py_k - i))

        return from_python(result)


class RandomPrime(Builtin):
    """
    <dl>
    <dt>'RandomPrime[{$imin$, $imax}]'
        <dd>gives a random prime between $imin$ and $imax$.
    <dt>'RandomPrime[$imax$]'
        <dd>gives a random prime between 2 and $imax$.
    <dt>'RandomPrime[$range$, $n$]'
        <dd>gives a list of $n$ random primes in $range$.
    </dl>

    >> RandomPrime[{14, 17}]
     = 17

    >> RandomPrime[{14, 16}, 1]
     : There are no primes in the specified interval.
     = RandomPrime[{14, 16}, 1]

    >> RandomPrime[{8,12}, 3]
     = {11, 11, 11}

    >> RandomPrime[{10,30}, {2,5}]
     = ...

    #> RandomPrime[{10,12}, {2,2}]
     = {{11, 11}, {11, 11}}

    #> RandomPrime[2, {3,2}]
     = {{2, 2}, {2, 2}, {2, 2}}
    """

    messages = {
        'posdim': (
            'The dimensions parameter `1` is expected to be a positive '
            'integer or a list of positive integers.'),
        'noprime': 'There are no primes in the specified interval.',
        'prmrng': ('First argument `1` is not a positive integer or a list '
                   'of two positive integers.'),
        'posint': ('The paramater `1` describing the interval is expected to '
                   'be a positive integer.'),
    }

    rules = {
        'RandomPrime[imax_?NotListQ]': 'RandomPrime[{1, imax}, 1]',
        'RandomPrime[int_?ListQ]': 'RandomPrime[int, 1]',
        'RandomPrime[imax_?ListQ, n_?ArrayQ]': (
            'ConstantArray[RandomPrime[imax, 1], n]'),
        'RandomPrime[imax_?NotListQ, n_?ArrayQ]': (
            'ConstantArray[RandomPrime[{1, imax}, 1], n]'),
    }

    # TODO: Use random state as in other randomised methods within mathics

    def apply(self, interval, n, evaluation):
        'RandomPrime[interval_?ListQ, n_]'

        if not isinstance(n, Integer):
            evaluation.message('RandomPrime', 'posdim', n)
            return
        py_n = n.to_python()

        py_int = interval.to_python()
        if not (isinstance(py_int, list) and len(py_int) == 2):
            evaluation.message('RandomPrime', 'prmrng', interval)

        imin, imax = min(py_int), max(py_int)
        if imin <= 0 or not isinstance(imin, int):
            evaluation.message('RandomPrime', 'posint', interval.leaves[0])
            return

        if imax <= 0 or not isinstance(imax, int):
            evaluation.message('RandomPrime', 'posint', interval.leaves[1])
            return

        try:
            if py_n == 1:
                return from_python(sympy.ntheory.randprime(imin, imax + 1))
            return from_python([sympy.ntheory.randprime(imin, imax + 1)
                                for i in range(py_n)])
        except ValueError:
            evaluation.message('RandomPrime', 'noprime')
            return


class Quotient(Builtin):
    '''
    <dl>
    <dt>'Quotient[m, n]'
      <dd>computes the integer quotient of $m$ and $n$.
    </dl>

    >> Quotient[23, 7]
     = 3

    #> Quotient[13, 0]
     : Infinite expression Quotient[13, 0] encountered.
     = ComplexInfinity
    #> Quotient[-17, 7]
     = -3
    #> Quotient[-17, -4]
     = 4
    #> Quotient[19, -4]
     = -5
    '''

    attributes = ('Listable', 'NumericFunction')

    messages = {
        'infy': 'Infinite expression `1` encountered.',
    }

    def apply(self, m, n, evaluation):
        'Quotient[m_Integer, n_Integer]'
        py_m = m.get_int_value()
        py_n = n.get_int_value()
        if py_n == 0:
            evaluation.message('Quotient', 'infy', Expression('Quotient', m, n))
            return Symbol('ComplexInfinity')
        return Integer(py_m // py_n)
