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
        ['T', '|', 'R'],
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
if __name__ == "__main__":
    # generate
    n = 1000
    cfg = CFG(regex_rules, 'R')
    output_file_path = 'data/regex.txt'
    output = []
    for i in range(n):
        regex = cfg.sample(max_depth=15)
        output.append([i+1, regex])

    with open(output_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for row in output:
            writer.writerow(row)

    # test
    input_file_path = 'data/regex.txt'
    output_file_path = 'results/regex_results.txt'
    output = []
    with open(input_file_path, 'r', newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        for row in csv_reader:
            # print(row[0], row[1])
            regex_str = row[1]
            nfa = NFA.from_regex(regex_str)
            my_dfa = DFA.from_nfa(nfa, minify=True)

            membership = d.decide_CRASP_membership(my_dfa)
            # print(f"CRASP membership: {membership}")
            output.append([row[0], row[1], membership])

    with open(output_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        for row in output:
            writer.writerow(row)