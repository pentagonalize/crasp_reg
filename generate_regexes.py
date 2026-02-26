import pysemigroup
import itertools
from graphviz import Source
import csv
import dfa
from dfa import generate_dfa_diagram
from automata.fa.dfa import DFA
from automata.fa.nfa import NFA
import networkx as nx
from sympy import Matrix
import matplotlib.pyplot as plt
import decider as d
import argparse
import random
from typing import Dict, List
import csv
import os

class CFG:
    """Context-Free Grammar sampler"""

    def __init__(self, rules: Dict[str, List[List[str]]], start_symbol: str = 'S'):
        self.rules = rules
        self.start_symbol = start_symbol

    def sample(self, symbol: str = None, max_depth: int = 50) -> str:
        if symbol is None:
            symbol = self.start_symbol

        if symbol not in self.rules:
            return symbol

        if max_depth <= 0:
            # When depth exhausted, force terminal-only productions
            terminal_prods = [p for p in self.rules[symbol]
                            if all(s not in self.rules for s in p)]
            if terminal_prods:
                production = random.choice(terminal_prods)
            else:
                # If no pure terminal production, pick shortest
                production = min(self.rules[symbol], key=len)

            # Expand remaining symbols with depth=0 (will force terminals)
            result = []
            for sym in production:
                if sym not in self.rules:
                    result.append(sym)
                else:
                    result.append(self.sample(sym, 0))
            return ''.join(result)

        production = random.choice(self.rules[symbol])

        result = []
        for sym in production:
            expanded = self.sample(sym, max_depth - 1)
            result.append(expanded)

        return ''.join(result)

    def sample_n(self, n: int, symbol: str = None) -> List[str]:
        return [self.sample(symbol) for _ in range(n)]


regex_rules = {
    # Regex with alternation
    'R': [
        ['T'],
        ['T', '+', 'R'],
    ],

    # Terms (concatenation)
    'T': [
        ['F'],
        ['F', 'T'],
    ],

    # Factors (atoms with optional star)
    'F': [
        ['A'],
        ['A', '*'],
    ],

    # Atoms (base symbols or grouped expressions)
    'A': [
        ['a'],
        ['b'],
        ['c'],
        ['d'],
        ['e'],
        ['(', 'R', ')'],
    ],
}


def check_R_infinity(regex, alphabet):
    # Create a semigroup from the regex
    language = pysemigroup.RegularLanguage.from_easy_regex(regex, alphabet)

    # get the syntactic monoid of the language
    semigroup = language.syntactic_monoid()

    # dot = semigroup.graphviz_string()
    # # replace print(dot) with:
    # Source(dot).view()  # writes a temp file and opens the default viewer

    # print("Elements of the semigroup:", list(semigroup.elements()))
    # print("Idempotents of the semigroup:", list(semigroup.idempotents()))

    two_idempotents_in_the_same_R_class = False
    for el1, el2 in itertools.combinations(list(semigroup.idempotents()), 2):
        if semigroup.R_class_of_element(el1) == semigroup.R_class_of_element(el2):
            two_idempotents_in_the_same_R_class = True
            # print(f"Two idempotents in the same R-class: {el1}, {el2}")
            break
    return not two_idempotents_in_the_same_R_class

if __name__ == "__main__":

    # generate
    n = 1000
    cfg = CFG(regex_rules, 'R')
    os.makedirs('data', exist_ok=True)
    os.makedirs('results', exist_ok=True)

    data_path = 'data/regex.txt'
    output_rows = []
    for i in range(n):
        regex = cfg.sample(max_depth=15)
        output_rows.append([i + 1, regex])

    with open(data_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for row in output_rows:
            writer.writerow(row)

    # classify into three bins:
    #  - C-RASP
    #  - R_infinity but not C-RASP
    #  - not in R_infinity
    input_file_path = data_path
    crasp_path = 'results/crasp.txt'
    r_infty_not_crasp_path = 'results/r_infinity_not_crasp.txt'
    not_r_infty_path = 'results/not_r_infinity.txt'
    errors_path = 'results/classify_errors.txt'

    # alphabet used by the regex generator
    alphabet = set(['a', 'b', 'c', 'd', 'e'])

    counts = {'crasp': 0, 'r_infty_not_crasp': 0, 'not_r_infty': 0, 'errors': 0}

    with open(input_file_path, 'r', newline='') as csvfile, \
         open(crasp_path, 'w', newline='') as crasp_file, \
         open(r_infty_not_crasp_path, 'w', newline='') as rinf_file, \
         open(not_r_infty_path, 'w', newline='') as notr_file, \
         open(errors_path, 'w', newline='') as err_file:

        reader = csv.reader(csvfile)
        crasp_writer = csv.writer(crasp_file)
        rinf_writer = csv.writer(rinf_file)
        notr_writer = csv.writer(notr_file)
        err_writer = csv.writer(err_file)

        # header rows (optional)
        crasp_writer.writerow(["index", "regex", "C-RASP", "R_infinity"])
        rinf_writer.writerow(["index", "regex", "C-RASP", "R_infinity"])
        notr_writer.writerow(["index", "regex", "C-RASP", "R_infinity"])
        err_writer.writerow(["index", "regex", "error"])

        for row in reader:
            try:
                idx = row[0]
                regex_str = row[1]

                # build NFA/DFA
                nfa = NFA.from_regex(regex_str)
                my_dfa = DFA.from_nfa(nfa, minify=True)

                # check C-RASP membership
                try:
                    membership = d.decide_CRASP_membership(my_dfa)
                except Exception as e:
                    # If decider fails, record error and skip
                    err_writer.writerow([idx, f"C-RASP decider error: {e}"])
                    counts['errors'] += 1
                    continue

                # check R_infinity
                try:
                    r_infty = check_R_infinity(regex_str, alphabet)
                except Exception as e:
                    err_writer.writerow([idx, f"R_infinity check error: {e}"])
                    counts['errors'] += 1
                    continue

                # binning
                if membership:
                    crasp_writer.writerow([idx, regex_str, True, r_infty])
                    counts['crasp'] += 1
                elif r_infty:
                    rinf_writer.writerow([idx, regex_str, False, True])
                    counts['r_infty_not_crasp'] += 1
                else:
                    notr_writer.writerow([idx, regex_str, False, False])
                    counts['not_r_infty'] += 1

            except Exception as e:
                # catch any unexpected parsing/IO error for the row
                err_writer.writerow(row + [f"unexpected error: {e}"])
                counts['errors'] += 1

    # write a small summary
    summary_path = 'results/classify_summary.txt'
    with open(summary_path, 'w') as sf:
        sf.write(f"Total generated: {n}\n")
        sf.write(f"C-RASP: {counts['crasp']}\n")
        sf.write(f"R_infinity but not C-RASP: {counts['r_infty_not_crasp']}\n")
        sf.write(f"Not R_infinity: {counts['not_r_infty']}\n")
        sf.write(f"Errors: {counts['errors']}\n")