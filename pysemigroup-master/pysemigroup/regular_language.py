#-*- coding:utf-8 -*-
import re
from .automata import Automaton
from .transition_semigroup import TransitionSemiGroup
_star = "_star"
r"""
%skip
EXAMPLES::

    sage: from pysemigroup import RegularLanguage
    sage: L = RegularLanguage("((a+b+c)^3)^_star*a*c*((a+c)^4)^_star*a*c*b*a*((a+b+c)^5)^_star",["a","b","c"])
    sage: L
    Regular language: ((a+b+c)^3)*ac((a+c)^4)*acba((a+b+c)^5)* over alphabet set(['a', 'c', 'b'])
    sage: A = L.automaton()
    sage: A
    Automaton of 76 states
    sage: B = A.minimal_automaton()                           #not tested - 71s
    sage: C = A.minimal_automaton(algorithm="Brzozowski")     #not tested - 63s
    sage: D = A.minimal_automaton(algorithm="Moore")          #not tested - 54s
    sage: B                         #not tested                                         
    Automaton of 2402 states         #not tested
    sage: C         #not tested
    Automaton of 2425 states         #not tested
    sage: D         #not tested
    Automaton of 2425 states         #not tested
    sage: A.is_accepted('aaaacaaaaacbaaaaaa')         #not tested
    True         #not tested
    sage: B.is_accepted('aaaacaaaaacbaaaaaa')         #not tested
    False         #not tested
    sage: C.is_accepted('aaaacaaaaacbaaaaaa')         #not tested
    True         #not tested
    sage: D.is_accepted('aaaacaaaaacbaaaaaa')         #not tested
    True

"""
from .automata import *
import itertools                    
    
class RegularLanguage:    
    def __init__(self, regex, letters=None):
        r"""
        INPUT:

        - ``regex`` - string different from the empty word
        - ``alphabet`` - set of atomic letters
        """
        self._regex = regex
        if letters is None:
            s = regex
            s = s.replace("(","")  
            s = s.replace(")","")        
            s = s.replace("^","")  
            s = s.replace("x","")        
            s = s.replace("*","")  
            s = s.replace("+","")
            s = s.replace(" ","")
            for x in range(10):
                s = s.replace(str(x),"")

            letters = set()
            for x in s:
                letters.add(x)
        self._letters = set(letters)

    @classmethod
    def from_easy_regex(cls, regex,A=None):
            r"""
            Generate a regular language using a simplify regex syntax. 
            The alphabet is automatically compute. To use with caution, 
            since in a lot of cases, it will not produce the desire regular language.
            INPUT:

            - ``alphabet`` -- list of letters 

            OUTPUT:

                Regular Languages over automatically compute alphabet

            EXAMPLES::

                sage: from pysemigroup import RegularLanguage
                sage: L = RegularLanguage.from_easy_regex("(abba)*")
                sage: RegularLanguage.from_easy_regex("(abc)*+(ad)*")
                Regular language: (abc)*+(ad)* over alphabet set(['a', 'c', 'b', 'd'])
                sage: RegularLanguage.from_easy_regex("(ab)^3+(cd)*")
                Regular language: (ab)^3+(cd)* over alphabet set(['a', 'c', 'b', 'd'])
            
                sage: RegularLanguage.from_easy_regex("(xx)*+(xd)*")
                Regular language: (xx)*+(xd)* over alphabet set(['x', 'd'])

            """
            if not(A):
                s = str(regex)
                s = s.replace("(","")  
                s = s.replace(")","")        
                s = s.replace("^","")  
                s = s.replace("_star","")        
                s = s.replace("*","")  
                s = s.replace("+","")
                s = s.replace("-","")
                s = s.replace(" ","")
                s = s.replace("A","")        
                for x in range(10):
                    while (str(x) in s):           
                        s = s.replace(str(x),"")
                A = set(s)
            alph_str = "("
            for a in A:
                alph_str = alph_str + a + "+"
            alph_str = alph_str[0:len(alph_str)-1] + ")"

            regex = regex.replace("A",alph_str)
            regex = regex.replace("1","(_empty_word)")  
            for x in A:
               for y in A:
                   regex = regex.replace(x+y,x+"."+y)
            regex = regex.replace("*","**_star")  
            regex = regex.replace(")(",")*(")
            regex = regex.replace("_star(","_star*(")
            for x in A:
               regex = regex.replace(x+"(",x+"*(")
               regex = regex.replace(")"+x,")*"+x)
               regex = regex.replace("_star"+x,"_star*"+x)
               regex = regex.replace(x+x,x+"*"+x)
            regex = regex.replace(".","*")           
            return RegularLanguage(regex,letters=A)
    def letters(self):
        r"""
        """
        return set(self._letters)
        
    def __repr__(self):
        r"""
        String representation

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: RegularLanguage("a*b", ['a', 'b'])
            Regular language: ab over alphabet set(['a', 'b'])

        """
        s = self._regex
        s = s.replace("*","")
        s = s.replace("_star","*")
        s = s.replace("1e",'1')
        return "Regular language: %s over alphabet %s" % (s, self._letters)
        
    def __iter__(self):
        r"""
        Return an iterator of the language.
        
        NOTE::
         
            Exponential complexity. TODO: improve it.
        
        EXAMPLES::
        
            sage: from pysemigroup import RegularLanguage
            sage: L = RegularLanguage("(a*b)^_star",["a","b"])
            sage: it = iter(L)
            sage: for _ in range(6): next(it)
            ''
            'ab'
            'abab'
            'ababab'
            'abababab'
            'ababababab'
        """
        letters = self._letters
        for n in itertools.count():
            for p in itertools.product(*(letters,)*n):
                if p in self:
                    yield ''.join(p).replace(" ","")
        
    def __contains__(self, word):
        r"""
        Return whether word is in the language.
        
        EXAMPLES::
        
            sage: from pysemigroup import RegularLanguage
            sage: L = RegularLanguage("(a*b)^_star",["a","b"])
            sage: "ababababab" in L
            True
            sage: "abababaabab" in L
            False
            
        Other notation::
        
            sage: ("a","b") in L
            True
            sage: ["b","a"] in L
            False

        """
        return self.automaton().is_accepted(word)
    def __eq__(self, other):
        r"""
        Return wheter self is equal to other.

        INPUT:

        -  ``self`` -  regular language
        -  ``other`` -  regular language

        OUTPUT:

        boolean

        EXAMPLES::


            sage: from pysemigroup import RegularLanguage
            sage: R = RegularLanguage("a*b", ['a', 'b'])
            sage: S = RegularLanguage("b*a*a", ['a', 'b'])
            sage: R == S
            False
            
        ::
            
            sage: L = RegularLanguage("a+b+c",["a","b","c"])
            sage: L == L
            True
            sage: L2 = RegularLanguage("a+a+b+c",["a","b","c"])
            sage: L == L2
            True
            sage: L3 = RegularLanguage("a*a+a+b+c",["a","b","c"])
            sage: L == L3
            False

        """    
        R = (self-other)+(other-self)
        return R.is_empty()
        
    def is_empty(self):
        r"""
        Return wheter self is empty.

        INPUT:

        -  ``self`` -  regular language

        OUTPUT:

        boolean

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: R = RegularLanguage("a*b", ['a', 'b'])
            sage: R.is_empty()
            False
            sage: R2 = RegularLanguage("-((a+b)^_star)", ['a', 'b'])
            sage: R2.is_empty()
            True


        """    
        A = self.automaton()
        return not A.is_finite_state_reachable()

    def __mul__(self, other):
        r"""
        Return the concatenation of regular languages.

        INPUT:

        -  ``self`` -  regular language
        -  ``other`` -  regular language

        OUTPUT:

        RegularLanguage

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: R = RegularLanguage("a*b", ['a', 'b'])
            sage: S = RegularLanguage("b*a*a", ['a', 'b'])
            sage: R * S
            Regular language: (abbaa) over alphabet set(['a', 'b'])

        """
        regex = "(%s*%s)" % (self._regex, other._regex)
        letters = self.letters() | other.letters()
        return RegularLanguage(regex, letters=letters)

    def __add__(self, other):
        r"""
        Return the union of regular languages.

        INPUT:

        -  ``self`` -  regular language
        -  ``other`` -  regular language

        OUTPUT:

        RegularLanguage

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: R = RegularLanguage("a*b", ['a', 'b'])
            sage: S = RegularLanguage("b*a*a", ['a', 'b'])
            sage: R + S
            Regular language: (ab+baa) over alphabet set(['a', 'b'])

        """
        regex = "(%s+%s)" % (self._regex, other._regex)
        letters = self.letters() | other.letters()
        return RegularLanguage(regex, letters=letters)
    
    def __sub__(self, other):
        r"""
        Difference symetric of self by other language.

        OUTPUT:

        RegularLanguage

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: L = RegularLanguage("(a*b)^_star", ['a', 'b']) 
            sage: R = RegularLanguage("a*b", ['a', 'b'])     
            sage: L-R
            Regular language: (((ab)*)-(ab)) over alphabet set(['a', 'b'])

        """
        regex = "((%s)-(%s))" % (self._regex, other._regex)
        letters = self.letters()
        return RegularLanguage(regex, letters)

    def __neg__(self):
        r"""
        Complement of self             

        OUTPUT:

        RegularLanguage

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: L = RegularLanguage("(a*b)^_star", ['a', 'b'])
            sage: -L
            Regular language: -(ab)* over alphabet set(['a', 'b'])

        

        """
        regex = "-%s" % self._regex
        letters = self.letters()
        return RegularLanguage(regex, letters)

    def intersection(self, other):
        r"""
        Return the intersection of regular languages.

        INPUT:

        -  ``self`` -  regular language
        -  ``other`` -  regular language

        OUTPUT:

        RegularLanguage

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: L = RegularLanguage("(a*b)^_star", ['a', 'b']) 
            sage: R = RegularLanguage("a*b", ['a', 'b'])     
            sage: L.intersection(R)
            Regular language: (((ab)*)-(-(ab))) over alphabet set(['a', 'b'])
            

        """
        regex = "((%s)-(-(%s)))" % (self._regex, other._regex)
        letters = self.letters() | other.letters()
        return RegularLanguage(regex, letters)        

    def kleene_star(self):
        r"""
        Return the Kleen star of a regular language.

        OUTPUT:

        RegularLanguage

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: R = RegularLanguage("a*b", ['a', 'b'])
            sage: R.kleene_star()
            Regular language: (ab)* over alphabet set(['a', 'b'])

        """
        regex = "(%s)^_star" % self._regex
        letters = self.letters()
        return RegularLanguage(regex, letters)

    def __pow__(self, exponent):
        r"""
        Return the Kleen star or the integer power of a regular language.

        INPUT:

        - ``exponent`` -- integer or variable x

        OUTPUT:

        RegularLanguage

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage,_star
            sage: R = RegularLanguage("a*b", ['a', 'b'])
            sage: R^_star
            Regular language: (ab)* over alphabet set(['a', 'b'])
            sage: R^3
            Regular language: (ab)^3 over alphabet set(['a', 'b'])

        """
        regex = "(%s)^%s" % (self._regex, exponent)
        letters = self.letters()
        return RegularLanguage(regex, letters)

    def automaton(self):
        r"""
        Return an automaton of the regular language.

        OUTPUT :

        automaton

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: RegularLanguage("a*a",["a"]).automaton()
            Automaton of 4 states
            sage: RegularLanguage("a^_star",["a"]).automaton()
            Automaton of 2 states
            sage: RegularLanguage("a*b",["a","b"]).automaton()
            Automaton of 4 states
            sage: RegularLanguage("a+b",["a","b"]).automaton()
            Automaton of 4 states
            sage: RegularLanguage("a^5",["a"]).automaton()
            Automaton of 10 states
            sage: RegularLanguage("((a^5)+b*b)^_star",["a","b"]).automaton()
            Automaton of 14 states

        """
        if self._regex == "":
            return Automaton.from_empty_string(self.letters())        
        else:  
            for letter in self.letters():
                exec(letter+" = Automaton.from_letter(letter,alphabet=self.letters())")
            _empty_word = Automaton.from_empty_string(self.letters())
            _star = "_star"
            return eval(self._regex)

    def automaton_deterministic(self):
        r"""
        Return a deterministic automaton of the regular language.

        OUTPUT :

        automaton

        EXAMPLES::

            sage: from pysemigroup import RegularLanguage
            sage: RegularLanguage("a*a",["a"]).automaton_deterministic()
            Automaton of 4 states
            sage: RegularLanguage("a^_star",["a"]).automaton_deterministic()
            Automaton of 3 states
            sage: RegularLanguage("a*b",["a","b"]).automaton_deterministic()
            Automaton of 4 states
            sage: RegularLanguage("a+b",["a","b"]).automaton_deterministic()
            Automaton of 4 states
            sage: RegularLanguage("a^5",["a"]).automaton_deterministic()
            Automaton of 7 states
            sage: RegularLanguage("((a^5)+b*b)^_star",["a","b"]).automaton_deterministic()
            Automaton of 9 states

        """
        return self.automaton().deterministic_automaton(rename_states=True)

    def automaton_minimal_deterministic(self, algorithm=None):
        r"""
        Return the minimal deterministic automaton of the regular language.
        INPUT:

        - ``algorithm`` -- None, or "Brzozowski" or "Moore"

        OUTPUT:

        automaton

        EXAMPLES::
        
            sage: from pysemigroup import RegularLanguage
            sage: RegularLanguage("a*a",["a"]).automaton_minimal_deterministic()
            Automaton of 4 states
            sage: RegularLanguage("a*b",["a","b"]).automaton_minimal_deterministic()
            Automaton of 4 states
            sage: RegularLanguage("a+b",["a","b"]).automaton_minimal_deterministic()
            Automaton of 3 states
            sage: RegularLanguage("a^5",["a"]).automaton_minimal_deterministic()
            Automaton of 7 states
            sage: RegularLanguage("((a^5)+b*b)^_star",["a","b"]).automaton_minimal_deterministic()
            Automaton of 7 states 
            sage: RegularLanguage("a^_star",["a"]).automaton_minimal_deterministic()
            Automaton of 1 states

        """
        return self.automaton().minimal_automaton(algorithm=algorithm)
  
    def syntactic_semigroup(self):
        return TransitionSemiGroup(self.automaton_minimal_deterministic(),monoid=False)

    def syntactic_monoid(self):
        return TransitionSemiGroup(self.automaton_minimal_deterministic())
  
    def view(self):
        self.syntactic_monoid().view()    
        self.automaton_minimal_deterministic().view()    

    def is_equation_satisfied(self, eq, variables, monoid=True, verbose=False):
        r"""
        EXAMPLES:

           

        """
        if monoid:
            return self.syntactic_monoid().is_equation_satisfied(eq,variables,verbose=verbose)
        else:                                        
            return self.syntactic_semigroup().is_equation_satisfied(eq,variables,verbose=verbose)
    
                   
