from .automata import Automaton
def prefix_set(words):
    prefix = set()
    for x in words:
        prefix.update([x[:i] for i in range(len(x)+1)])
    return prefix

def prefix_automaton(words, alphabet=None):
    if len(words) == 0:
        raise ValueError("Set of words should not be empty")
    if alphabet is None:
        alphabet = set()
        for word in words:
            alphabet.update(word)
    prefix = prefix_set(words)
    transitions = {}
    init = min(prefix, key=lambda e:len(e))
    for word in prefix:
        for letter in alphabet:
            if word+letter in prefix:
                transitions[(frozenset([word]), letter)] = [frozenset([word + letter])]
    return Automaton(transitions, initial_states=[frozenset([init])], final_states=[frozenset([word]) for word in words],  alphabet=alphabet)

def rpni_merge(automaton, state1, state2):   
    transitions = {}
    old_transitions = automaton._transitions
    old_states = automaton._states
    for ostate in old_states:
        if len(state1.intersection(ostate))>0 :
            bag1 = ostate
        if len(state2.intersection(ostate))>0:
            bag2 = ostate
    merge_state = bag1.union(bag2)
    states = set(old_states)
    states.remove(bag1)
    if bag2 in states:
        states.remove(bag2)
    states.add(merge_state)
    for transition in old_transitions:
        target = old_transitions[transition][0]
        source = transition[0]
        letter = transition[1]
        if len(merge_state.intersection(target)) > 0:
            target = merge_state
        if len(merge_state.intersection(transition[0]))>0:
            source = merge_state
        if (source, letter) not in transitions:
            transitions[(source, letter)] = []
        transitions[(source, letter)].append(target)
    initial_states = [state for state in states
                      if len(list(filter(lambda e: len(e.intersection(state))>0, automaton._initial_states ))) > 0]
    final_states = [state for state in states
                      if len(list(filter(lambda e: len(e.intersection(state))>0, automaton._final_states))) > 0]    
    return Automaton(transitions, initial_states=initial_states, final_states=final_states)

def rpni_simplify(automaton):
    transitions = automaton._transitions
    A = automaton
    while not A.is_deterministic():
        states_to_merge = list(filter(lambda e: len(e)>1, A._transitions.values())).pop()
        A = rpni_merge(A, states_to_merge[0], states_to_merge[1])        
    return A

def rpni(positive_words, negative_words, alphabet=None, rename=True):
    positif_aut = prefix_automaton(positive_words, alphabet=alphabet)
    negative_aut = prefix_automaton(negative_words, alphabet=alphabet).minimal_automaton()
    ordered_states = sorted(sorted(positif_aut._states, key=lambda a:list(a)[0]), key=lambda a: len(list(a)[0]))
    for state1 in ordered_states:
        for state2 in ordered_states[:ordered_states.index(state1)]:
            new_aut = rpni_simplify(rpni_merge(positif_aut, state1, state2))            
            if not (new_aut.intersection(negative_aut).is_finite_state_reachable()):
                positif_aut = new_aut
    if rename: 
        positif_aut.rename_states()
    return positif_aut
