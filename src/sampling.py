import itertools
import random
import logging
from pathlib import Path
from collections import defaultdict
from typing import Tuple, Hashable, Set, Dict, List
import decimal as d
import pprint

import networkx as nx

PauliTerm = str  # {I, X, Y, Z}^n or {X, Y, Z}^n
PauliOp = str    # { I, X, Y, Z} or {X, Y, Z}

PAULI_COLORS = {
    'I': 'gray',
    'X': 'red',
    'Y': 'darkgreen',
    'Z': 'blue'
}

LOG_LEVELS = {
    0: logging.WARN,
    1: logging.INFO,
    2: logging.DEBUG
}


def write_to_dot(graph: nx.MultiDiGraph, filename: str) -> None:
    """There is a bug in GraphViz 2.40 that ignores some edges in multi-graphs
    so we handle the to-dot functionality ourselves"""
    lines = [
        'digraph "" {',
        'node [shape=circle, label="\\N"];',
        'edge [fontsize=10];',
        '-1 [label=1,shape=square];',
    ]
    for u, v, data in graph.edges(data=True):
        lines.append('{u} -> {v} [color={color}, fontcolor={color}, label="{weight}{pauli}{virtual}({paths})"];'.format(
            u=u, v=v,
            color=PAULI_COLORS[data['pauli']],
            weight=(f'{data["weight"]:.8f}' if data['weight'] != 1 else ''),
            pauli=data['pauli'],
            virtual=('^' if data.get('virtual', False) else ''),
            paths=data.get('paths', '-')
        ))

    lines.append("}")
    with Path(filename).open('w') as f:
        logging.debug(f'writing graph to {filename}')
        f.write('\n'.join(lines))


def write_hamiltonian(H: Dict[PauliTerm, d.Decimal], filename: str) -> None:
    lines = []
    for p, c in H.items():
        lines.append(p)
        lines.append(f'{c}')

    with Path(filename).open('w') as f:
        f.write('\n'.join(lines))


def has_pauli_out_edge(graph: nx.MultiDiGraph, node: Hashable, pauli: PauliTerm) -> bool:
    out_edges = graph.out_edges(node, data=True)
    for edge in out_edges:
        if edge[2]['pauli'] == pauli:
            return True
    return False


def get_pauli_weight(graph: nx.MultiDiGraph, node: Hashable, pauli: PauliTerm) -> Tuple[Hashable, d.Decimal]:
    out_edges = graph.out_edges(node, data=True)
    for edge in out_edges:
        if edge[2]['pauli'] == pauli:
            return edge[1], edge[2]['weight']
    raise ValueError('Node does not have a pauli-{} edge'.format(pauli))


def compare_nodes(graph: nx.MultiDiGraph, node_a: Hashable, node_b: Hashable) -> bool:
    if graph.nodes[node_a]['level'] != graph.nodes[node_b]['level']:
        return False
    if graph.nodes[node_a]['hash'] != graph.nodes[node_b]['hash']:
        return False

    for pauli in "IXYZ":
        has_pauli_a = has_pauli_out_edge(graph, node_a, pauli)
        has_pauli_b = has_pauli_out_edge(graph, node_b, pauli)
        if has_pauli_a != has_pauli_b:
            return False
        if has_pauli_a:
            target_a, weight_a = get_pauli_weight(graph, node_a, pauli)
            target_b, weight_b = get_pauli_weight(graph, node_b, pauli)
            if target_a != target_b or weight_a != weight_b:
                return False
    return True


def count_paths(graph: nx.MultiDiGraph, curr_node: Hashable, counter_dict: Dict[Hashable, int]) -> int:
    """
    count the number of paths passing thru a node
    """
    if curr_node in counter_dict:
        return counter_dict[curr_node]

    out_edges = graph.out_edges(curr_node, data=True)
    if out_edges:
        _count = 0
        for u, v, data in out_edges:
            if v not in counter_dict:
                counter_dict[v] = count_paths(graph, v, counter_dict)
            _count += counter_dict[v]
            data["paths"] = counter_dict[v]
        counter_dict[curr_node] = _count
        return _count
    else:  # no out edges
        counter_dict[curr_node] = 1
    return 1


def sample_from_paths_uniform(graph: nx.MultiDiGraph, source: Hashable = 0) -> Tuple[d.Decimal, PauliTerm]:
    current_node = source
    string = []
    value = d.Decimal(1.0)
    while current_node != -1:
        out_edges = list(graph.out_edges(current_node, data=True))
        u, v, data = random.choices(out_edges, weights=[dd['paths'] for _, _, dd in out_edges])[0]
        value *= data['weight']
        if data['pauli'] != 'I':
            string.append(data['pauli'])
        else:
            string.append(random.choice('XYZ'))
        current_node = v

    return value, ''.join(string)


def sample_from_nodes_uniform(graph: nx.MultiDiGraph, source: Hashable = 0) -> Tuple[d.Decimal, PauliTerm]:
    current_node = source
    string = []
    value = d.Decimal(1.0)
    while current_node != -1:
        out_edges = list(graph.out_edges(current_node, data=True))
        u, v, data = random.choice(out_edges)
        value *= data['weight']
        if data['pauli'] != 'I':
            string.append(data['pauli'])
        else:
            string.append(random.choice('XYZ'))
        current_node = v

    return value, ''.join(string)


def sample_from_coefficients(graph: nx.MultiDiGraph, source: Hashable = 0) -> Tuple[d.Decimal, PauliTerm]:
    current_node = source
    string = []
    value = d.Decimal(1.0)
    while current_node != -1:
        out_sum = sum(w for u, v, w in graph.out_edges(current_node, data='weight'))
        assert out_sum < d.Decimal(1.01), f'Failed {out_sum} < 1.01!'
        threshold = d.Decimal(random.random())
        sum_of_weight = d.Decimal(0)
        for u, v, data in graph.out_edges(current_node, data=True):
            sum_of_weight += data['weight']
            if sum_of_weight > threshold:
                value *= data['weight']
                if data['pauli'] != 'I':
                    string.append(data['pauli'])
                else:
                    string.append(random.choice('XYZ'))
                current_node = v
                break

    return value, ''.join(string)


def read_file_to_graph(input_file: Path) -> Tuple[nx.MultiDiGraph, Set]:
    logging.info("Reading file {}".format(input_file))
    hamiltonian = {}
    with input_file.open('r') as fp:
        for string, coefficient in itertools.zip_longest(fp, fp):
            string: str = string.upper().strip()
            coefficient = coefficient.strip()
            if string.startswith('END') or coefficient.startswith('END'):
                logging.info('found END marker in input file')
                break
            assert string not in hamiltonian, f'{string} was found a second time'
            hamiltonian[string] = d.Decimal(coefficient)
    graph, possible_paths = read_dict_to_graph(hamiltonian)
    graph.name = input_file.stem
    return graph, possible_paths


def read_file_to_graph_with_ldf_info(input_file: Path) -> Tuple[nx.MultiDiGraph, Set]:
    logging.info("Reading file {}".format(input_file))
    hamiltonian = defaultdict(int)
    current_ldf_term = None
    with input_file.open('r') as fp:
        for string, coefficient in itertools.zip_longest(fp, fp):
            string: str = string.upper().strip()
            coefficient: str = coefficient.strip()
            if string.startswith('END') or coefficient.startswith('END'):
                logging.info('found END marker in input file')
                break
            if all(p == 'I' for p in string):
                continue
            coefficient: d.Decimal = abs(d.Decimal(coefficient))
            if coefficient == 0:
                current_ldf_term = string
            hamiltonian[current_ldf_term] += d.Decimal(coefficient)
    pprint.pprint(hamiltonian)
    graph, possible_paths = read_dict_to_graph(dict(hamiltonian))
    graph.name = input_file.stem
    return graph, possible_paths


def reduce_hamiltonian_dict_uniformly(hamiltonian: Dict[PauliTerm, d.Decimal]) -> Dict[PauliTerm, d.Decimal]:
    def mp(p1: PauliTerm, p2: PauliTerm) -> PauliTerm:
        term: List[str] = []
        for c1, c2 in zip(p1, p2):
            if c1 == 'I':
                term.append(c2)
            elif c2 == 'I':
                term.append(c1)
            elif c1 == c2:
                term.append(c1)
            else:
                raise RuntimeError(f'Mismatch during merge: {c1} vs. {c2}.')
        return ''.join(term)

    keys: List[PauliTerm] = sorted(list(hamiltonian.keys()), reverse=False, key=lambda p: (-p.count('I'), p))
    while keys:
        pauli1 = keys.pop(0)
        coeff1 = hamiltonian[pauli1]

        mergers = []

        for pauli2, coeff2 in sorted(hamiltonian.items(), reverse=False, key=lambda e: e[0]):
            if pauli1 != pauli2 and all(c1 == c2 or c1 == 'I' or c2 == 'I' for c1, c2 in zip(pauli1, pauli2)):
                mergers.append(pauli2)
                break

        if mergers:
            hamiltonian.pop(pauli1)

            for pauli2 in mergers:
                coeff2 = hamiltonian.pop(pauli2)
                keys.remove(pauli2)
                new_pauli = mp(pauli1, pauli2)
                if new_pauli not in hamiltonian:
                    hamiltonian[new_pauli] = coeff2 + coeff1 / len(mergers)
                else:
                    hamiltonian[new_pauli] += coeff2 + coeff1 / len(mergers)
                if new_pauli not in keys:
                    keys.append(new_pauli)
        keys.sort(reverse=True, key=lambda p: (-p.count('I'), p))
    return hamiltonian


def reduce_hamiltonian_dict_grouping(hamiltonian: Dict[PauliTerm, d.Decimal]) -> Dict[PauliTerm, d.Decimal]:
    """Coefficients at this point are assumed to be positive"""
    def merge_paulis(p1: PauliTerm, p2: PauliTerm) -> PauliTerm:
        """Merge Pauli terms such that the resulting term covers both"""
        term: List[str] = []
        for c1, c2 in zip(p1, p2):
            if c1 == 'I':
                term.append(c2)
            elif c2 == 'I':
                term.append(c1)
            elif c1 == c2:
                term.append(c1)
            else:
                raise RuntimeError(f'Mismatch during merge: {c1} vs. {c2}.')
        return ''.join(term)

    converged = False

    while not converged:
        converged = True
        groups = defaultdict(list)
        for pauli1, pauli2 in itertools.product(hamiltonian.keys(), repeat=2):
            if pauli1 != pauli2 and all(c1 == c2 or c1 == 'I' or c2 == 'I' for c1, c2 in zip(pauli1, pauli2)):
                groups[pauli1].append(pauli2)

        if groups:
            pauli_orig, pauli_targets = max(groups.items(), key=lambda e: (len(e[1]), e[0].count('I')))

            coeff_orig = hamiltonian.pop(pauli_orig)

            for pauli_t in pauli_targets:
                coeff_t = hamiltonian.pop(pauli_t)
                new_pauli = merge_paulis(pauli_orig, pauli_t)
                if new_pauli in hamiltonian:
                    hamiltonian[new_pauli] += coeff_t + coeff_orig / len(pauli_targets)
                else:
                    hamiltonian[new_pauli] = coeff_t + coeff_orig / len(pauli_targets)

            converged = False

    return hamiltonian



def read_dict_to_graph(hamiltonian: Dict[PauliTerm, d.Decimal]) -> Tuple[nx.MultiDiGraph, Set]:
    hamiltonian = {k: abs(v) for k, v in hamiltonian.items() if any(p != 'I' for p in k) and abs(v) > 0}
    #hamiltonian = reduce_hamiltonian_dict_uniformly(hamiltonian)
    hamiltonian = reduce_hamiltonian_dict_grouping(hamiltonian)
    #write_hamiltonian(hamiltonian, 'output/h.txt')
    graph = nx.MultiDiGraph()
    possible_paths = set()
    graph.add_node(0, level=0)
    graph.add_node(-1, level=-1, hash=1337)
    i = 0
    ignored = 0
    for pauli_str, coefficient in hamiltonian.items():
        pauli_str: str = pauli_str.upper().strip()
        if pauli_str.startswith('END'):
            logging.info('found END marker in input file')
            break
        if d.Decimal(abs(coefficient)) == 0:
            logging.debug(f'term {pauli_str} has too small coefficient of {coefficient} (term ignored)')
            ignored += 1
            continue
        possible_paths.add(pauli_str)
        current_node = 0
        for level, char in enumerate(pauli_str[:-1], 1):
            for u, v, data in graph.out_edges(current_node, data=True):
                if data['pauli'] == char:
                    current_node = v
                    break
            else:  # executed when no break in loop
                i += 1
                graph.add_node(i, level=level)
                graph.add_edge(current_node, i,
                               pauli=char,
                               key=char,
                               weight=d.Decimal(1.0))
                current_node = i
        graph.add_edge(current_node, -1,
                       pauli=pauli_str[-1],
                       key=pauli_str[-1],
                       weight=abs(d.Decimal(coefficient)))
    #print('ignored pauli terms = ', ignored)
    return graph, possible_paths


def renormalize_node(graph: nx.MultiDiGraph, node: Hashable, min_edge_weight: d.Decimal = d.Decimal(0)) -> None:
    out_edges = list(graph.out_edges(node, data=True))
    weight_sum = sum(edge[2]['weight'] for edge in out_edges)
    assert 1 <= len(out_edges) <= 4, f'len(out_edges) == {len(out_edges)} for node {node}'

    re_renormalize = False
    for u, v, data in out_edges:
        weight_before = data['weight']
        data['weight'] /= weight_sum

        if data['weight'] < min_edge_weight:
            data['weight'] = min_edge_weight
            re_renormalize = True

    if re_renormalize:
        renormalize_node(graph, node)
    else:
        if abs(sum(data['weight'] for u, v, data in out_edges) - 1) >= 0.0001:
            for u, v, data in graph.out_edges(node, data=True):
                print(f'{u} -{data["pauli"]}-> {v}    {data["weight"]} ({data.get("virtual", False)})')
            raise RuntimeError('Sum of weights not 1')


def renormalize_nodes_recursively(graph: nx.MultiDiGraph, node: Hashable, min_edge_weight: d.Decimal = d.Decimal(0)) -> None:
    out_edges = list(graph.out_edges(node, data=True))
    in_edges = list(graph.in_edges(node, data=True))
    weight_out_sum = sum(edge[2]['weight'] for edge in out_edges)
    assert 1 <= len(out_edges) <= 4, f'len(out_edges) == {len(out_edges)} for node {node}'

    re_renormalize = False
    for u, v, data in out_edges:
        data['weight'] /= weight_out_sum

        if data['weight'] < min_edge_weight:
            data['weight'] = min_edge_weight
            re_renormalize = True

    if re_renormalize:
        weight_out_sum_new = sum(edge[2]['weight'] for edge in out_edges)
        for u, v, data in out_edges:
            data['weight'] /= weight_out_sum_new
        weight_out_sum *= weight_out_sum_new

    for u, v, data in in_edges:
        data['weight'] *= weight_out_sum
    for u, v, data in in_edges:
        renormalize_nodes_recursively(graph, u)


def hash_node(graph: nx.MultiDiGraph, node: Hashable) -> int:
    out_edges = graph.out_edges(node, data=True)
    hash_dict = {}
    for u, v, data in out_edges:
        hash_dict[data['pauli']] = (data['weight'], graph.nodes[v]['hash'])
    return hash(frozenset(hash_dict))


def normalize_graph(graph: nx.MultiDiGraph, possible_paths: Set[PauliTerm]) -> None:
    """
        Given a graph with this function calculates the appropriate edge weights and eliminates the identity edges.
    """
    write_to_dot(graph, 'output/{}_init.dot'.format(graph.name))
    # Set up edge weights for for graph
    node_queue = list(set([u for u, v in graph.in_edges(-1)]))
    while node_queue:
        current_node = node_queue.pop(0)
        out_edges = graph.out_edges(current_node, data=True)
        weight_sum: d.Decimal = sum(data['weight'] for u, v, data in out_edges)
        for u, v, data in out_edges:
            data['weight'] /= weight_sum

        for u, v, data in graph.in_edges(current_node, data=True):
            data['weight'] *= weight_sum
            if u not in node_queue:
                node_queue.append(u)

    assert check_paths_and_graph(graph, possible_paths)
    write_to_dot(graph, 'output/{}_normalized.dot'.format(graph.name))

    merge_equivalent_node(graph)
    write_to_dot(graph, 'output/{}_mergefirst.dot'.format(graph.name))
    remove_identities(graph)
    remove_spurious_virtual_edges(graph)
    #prune_negligible_edges(graph, d.Decimal("1e-7"))  # this biases the probability distribution
    merge_equivalent_node(graph)

    count_paths(graph, 0, dict())


def merge_equivalent_node(graph: nx.MultiDiGraph) -> None:
    # normalize graph by merging similar nodes and calculating hashes of nodes AGAIN
    node_queue = list(set([u for u, v in graph.in_edges(-1)]))
    unique_table = defaultdict(dict)
    i = 0
    while node_queue:
        i += 1
        current_node: Hashable = node_queue.pop(0)
        graph.nodes[current_node]['hash'] = hash_node(graph, current_node)

        if (graph.nodes[current_node]['hash'] in unique_table[graph.nodes[current_node]['level']] and
                compare_nodes(graph, current_node, unique_table[graph.nodes[current_node]['level']][graph.nodes[current_node]['hash']])):
            unique_node = unique_table[graph.nodes[current_node]['level']][graph.nodes[current_node]['hash']]
            if current_node != unique_node:
                in_edges = list(graph.in_edges(current_node, data=True))
                graph.remove_node(current_node)
                for u, v, data in in_edges:
                    graph.add_edge(u, unique_node, key=data['pauli'], **data)
                    logging.debug(f'Added edge {u} -> {unique_node} (was {u} -> {v})')
                    renormalize_node(graph, u)
                    node_queue.append(u)
                assert current_node != 0
                node_queue = [n for n in node_queue if n != current_node]
                logging.debug('... removed {} since similar to {}'.format(current_node, unique_node))

        else:
            unique_table[graph.nodes[current_node]['level']][graph.nodes[current_node]['hash']] = current_node
            for edge in graph.in_edges(current_node, data=True):
                if edge[0] not in node_queue:
                    assert 'hash' not in graph[edge[0]]
                    # logging.debug('{} not in queue. Appending...'.format(edge[0]))
                    node_queue.append(edge[0])


def remove_identities(graph: nx.MultiDiGraph) -> None:
    # eliminate the identity for easy cases first: Either u-I->v is the only out edge from u
    # or there is another edge u-P->v covering for the identity as well (P in XYZ)
    ONE_THIRD = d.Decimal(1) / d.Decimal(3)
    for node in graph.nodes:
        out_edges = list(graph.out_edges(node, data=True))
        has_identity = any(data['pauli'] == 'I' for u, v, data in out_edges)

        if has_identity and len(out_edges) == 1:
            # identity is the only out-going edge. Split into virtual edges
            u, v, data = out_edges[0]
            assert data['weight'] == 1
            graph.add_edge(u, v, key='X', pauli='X', weight=ONE_THIRD, virtual=True)
            graph.add_edge(u, v, key='Y', pauli='Y', weight=ONE_THIRD, virtual=True)
            graph.add_edge(u, v, key='Z', pauli='Z', weight=ONE_THIRD, virtual=True)
            graph.remove_edge(u, v, key='I')
        elif has_identity and len(out_edges) > 1:
            # check if ordinary edge and virtual edge share source and target
            identity_v: Hashable = [v for u, v, data in out_edges if data['pauli'] == 'I'][0]
            else_v: List[Hashable] = [v for u, v, data in out_edges if data['pauli'] != 'I']
            if identity_v in else_v:
                i_weight = [data['weight'] for u, v, data in out_edges if data['pauli'] == 'I'][0]
                neighor_edges = sum(1 for u, v, data in out_edges if data['pauli'] != 'I' and v == identity_v)
                graph.remove_edge(node, identity_v, key='I')
                for u, v, data in graph.out_edges(node, data=True):
                    if v == identity_v:
                        data['weight'] += i_weight / neighor_edges

    # eliminate the remaining identities by merging nodes
    write_to_dot(graph, f'output/{graph.name}_preIE.dot')
    node_queue: List[Hashable] = list(set([u for u, v in graph.in_edges(-1)]))
    i: int = 0
    while node_queue:
        current_node: Hashable = node_queue.pop(0)

        for u, v in graph.in_edges(current_node):
            if u not in node_queue:
                node_queue.append(u)

        out_edges = list(graph.out_edges(current_node, data=True))
        has_identity = any(data['pauli'] == 'I' for u, v, data in out_edges)

        if not has_identity:
            continue

        assert sum(1 for u, v, data in out_edges if data['pauli'] == 'I') == 1

        # so the current node at this point has an `u-I->v` edge and at least one XYZ-edge
        # which the I-edge being the only one to `v`
        identity_e = [(u, v, data) for u, v, data in out_edges if data['pauli'] == 'I'][0]
        else_es: List = [(u, v, data) for u, v, data in out_edges if data['pauli'] != 'I']

        # else_e = else_es[0]  # TODO: How to chose the best edge?
        # else_e = max(else_es, key=lambda e: e[2]['weight'])
        else_e = min(else_es, key=lambda e: e[2]['weight'])
        #print(f'Merging {identity_e[0]: 4}-{identity_e[2]["pauli"]}->{identity_e[1]: 4} {identity_e[2]["weight"]:.3f} into {else_e[0]: 4}-{else_e[2]["pauli"]}->{else_e[1]: 4} {else_e[2]["weight"]:.3f} out of {len(else_es)}')
        merge_into(graph, current_node, 'I', identity_e[1], else_e[1])
        else_e[2]['weight'] += identity_e[2]['weight']

        graph.remove_edge(current_node, identity_e[1], key='I')

        outg_sum = sum(data['weight'] for u, v, data in graph.out_edges(current_node, data=True))
        assert abs(outg_sum - 1) < 0.001, f"outg_sum == {outg_sum}"

        renormalize_node(graph, current_node)
        #write_to_dot(graph, f'output/{graph.name}_ie{i:04}.dot')
        i += 1
    remove_spurious_roots(graph)
    for n in graph.nodes:
        if n != -1:
            renormalize_node(graph, n, d.Decimal('0.1'))


def merge_into(graph: nx.MultiDiGraph, parent: Hashable, ps_pauli: PauliTerm, source: Hashable, target: Hashable, level: int = 0) -> None:
    """
    Merging a node `source` into another node `target` (with a common parent node  `parent`) and a given
    Pauli `ps_pauli` to uniquely identify the correct edge between `parent` and `source`.

    The out-edges of `source` are iterated and individually merged into `target`.
    Each edge `e` (source, dest, pauli, weight) falls into one of the following cases:

    1. `target` has *no* out-edge with the Pauli of `e`.
        A new edge from (target, dest, pauli, weight) is added and the weights of out-edges of `target` are adjusted.
    2. `target` has an out-edge with the Pauli of `e` *and* both edges point to the same successor node.
        The weights of out-edges of `target` are adjusted.
    3. Otherwise. (`source` and `target` have out-edges with the same Pauli pointing to different nodes)
        Find `v` such that is is the successor of `target` with Pauli `pauli`.
        Recursively call merge_into(`graph`, `source`, `pauli`, `dest`, `v`).
    """
    source_out = {data['pauli']: (v, data) for u, v, data in graph.out_edges(source, data=True)}
    target_out = {data['pauli']: (v, data) for u, v, data in graph.out_edges(target, data=True)}

    #print(f'[level={level}] Merging {source}:{ps_pauli} into {target}')

    for source_p, (source_v, source_data) in source_out.items():
        if source_p not in target_out:
            # target does not have an out-edge with the same pauli as the source edge
            #print(f'[level={level}] {source_p} not in {target_out}')
            graph.add_edge(target, source_v, key=source_data['pauli'], **source_data)
            graph.nodes[target]['merge_target'] = True
            renormalize_node(graph, target)

        elif source_out[source_p][0] == target_out[source_p][0]:
            # source edge and target edge point to same node (may be the terminal node)
            #print(f'[level={level}] {source}-{source_p}->{source_out[source_p][0]} ({source_out[source_p][1]["weight"]:.4f})  into {target}-{source_p}->{target_out[source_p][0]} ({target_out[source_p][1]["weight"]:.4f})')
            #for u, v, data in graph.out_edges(target, data=True):
            #    print(f'          {u}-{data["pauli"]}->{v} ({data["weight"]:.8f})')
            #target_out[source_p][1]['weight'] += source_out[source_p][1]['weight']
            graph.nodes[target]['merge_target'] = True
            target_out[source_p][1]['virtual'] = target_out[source_p][1].get('virtual', False) and source_out[source_p][1].get('virtual', False)
        else:
            # we have to look deeper to be able to merge
            #print(f'[level={level}] else')
            assert source_p in source_out and source_p in target_out
            merge_into(graph, source, source_p, source_v, target_out[source_p][0], level+1)
            target_out[source_p][1]['weight'] += source_out[source_p][1]['weight']
            graph.nodes[target]['merge_target'] = True
            target_out[source_p][1]['virtual'] = target_out[source_p][1].get('virtual', False) and source_out[source_p][1].get('virtual', False)


def remove_spurious_virtual_edges(graph: nx.MultiDiGraph) -> None:
    """Remove virtual edges that are not required after normalization"""
    for u in graph.nodes():
        out_edges = list(graph.out_edges(u, data=True))
        num_v_edges = sum(1 for u, v, data in out_edges if data.get('virtual', False))
        num_o_edges = sum(1 for u, v, data in out_edges if not data.get('virtual', False))

        # node does not have outgoing virtual edges. Nothing to do
        if num_v_edges == 0:
            continue

        if num_o_edges > 0 and num_v_edges > 0:
            # node has virtual edges and at least an ordinary edge
            # remove all virtual edges
            for u, v, data in out_edges:
                if data.get('virtual', False):
                    graph.remove_edge(u, v, key=data['pauli'])
        elif num_v_edges > 0:
            # all outgoing edges are virtual
            # just take one edge and remove the rest
            assert all(data.get('virtual', False) for u, v, data in out_edges)
            for u, v, data in out_edges[1:]:
                out_edges[0][2]['weight'] += data['weight']
                #print(f'Removed virtual edge {u} -{data["pauli"]}-> {v} {data["weight"]}')
                graph.remove_edge(u, v, key=data['pauli'])
        renormalize_node(graph, u)

    remove_spurious_roots(graph)


def prune_negligible_edges(graph: nx.MultiDiGraph, prune_below: d.Decimal) -> None:
    """Remove edges with a weight below `prune_below`. Makes probability distribution biased"""
    for u, v, data in list(graph.edges(data=True)):
        if data['weight'] < prune_below:
            graph.remove_edge(u, v, key=data['pauli'])
        renormalize_node(graph, u)
    remove_spurious_roots(graph)


def remove_spurious_roots(graph: nx.MultiDiGraph, root: Hashable = 0) -> None:
    before, after = 0, graph.number_of_nodes()
    while before != after:
        before = graph.number_of_nodes()
        graph.remove_nodes_from([node for node, in_degree in graph.in_degree() if in_degree == 0 and node != root])
        after = graph.number_of_nodes()


def dfs(graph: nx.MultiDiGraph, current_node: Hashable, path: str, paths: Set[str]) -> Set[str]:
    out_edges = graph.out_edges(current_node, data='pauli')

    if out_edges:
        for u, v, pauli in out_edges:
            new_path = path + pauli
            paths = dfs(graph, v, new_path, paths)
    else:
        paths.add(path)
    return paths


def graph_to_dict(graph: nx.MultiDiGraph) -> Dict[PauliTerm, d.Decimal]:
    paths = dfs(graph, 0, '', set())

    new_dict = {}
    for p in paths:
        new_dict[p] = compute_covered_prob(graph, 0, p)
    return new_dict


def check_paths_and_graph(graph: nx.MultiDiGraph, possible_paths: Set[PauliTerm]) -> bool:
    source = 0
    paths_in_graph = dfs(graph, source, '', set())

    if possible_paths == paths_in_graph:
        return True
    logging.error('|P_possible| = {}'.format(len(possible_paths)))
    logging.error('|P_graph|    = {}'.format(len(paths_in_graph)))
    logging.error('|P_p - P_g|  = {} {!r:.50}...'.format(len(possible_paths - paths_in_graph), possible_paths - paths_in_graph))
    logging.error('|P_g - P_p|  = {} {!r:.50}...'.format(len(paths_in_graph - possible_paths), paths_in_graph - possible_paths))
    logging.error('|P_g & P_p|  = {}'.format(len(paths_in_graph & possible_paths)))
    logging.error('|P_g ^ P_p|  = {} {!r:.50}...'.format(len(paths_in_graph ^ possible_paths), paths_in_graph ^ possible_paths))


def paths_match(path_without_i: PauliTerm, path_with_i: PauliTerm) -> bool:
    if len(path_without_i) != len(path_with_i):
        return False

    for cwo, cwi in zip(path_without_i, path_with_i):
        if cwi == 'I':
            continue
        elif cwi != cwo:
            return False
    return True


def check_paths_and_graph_with_wildcard(graph: nx.MultiDiGraph, possible_paths: Set[PauliTerm]) -> bool:
    paths_in_graph = dfs(graph, 0, '', set())

    ret = True
    for ppath in possible_paths:
        for gpath in paths_in_graph:
            if paths_match(gpath, ppath):
                break
        else:  # read if not break: return False
            logging.error('path {} not found in graph'.format(ppath))
            ret = False
    return ret


def print_some_stats(graph: nx.MultiDiGraph):
    level_count = defaultdict(int)
    for node, level in graph.nodes(data='level'):
        level_count[level] += 1
    for level, count in level_count.items():
        logging.debug('level {} has {} nodes'.format(level, count))


def compute_covered_prob(graph: nx.MultiDiGraph, curr_node: Hashable, pauliString: PauliTerm) -> d.Decimal:
    """
        (Newly added by Rudy Raymond)
        Compute the probability that the pauliString is covered by the walk on graph starting from curr_node
        The computation of the probability is similar to dfs
        ASSUMING: only pauli edges and no Identity edges in the graph
    """
    if len(pauliString) <= 0:
        return d.Decimal(1)

    out_edges = list(graph.out_edges(curr_node, data=True))

    if not out_edges:
        raise ValueError(f'No edges but pauli string "{pauliString}" not consumed')


    if pauliString[0] in 'XYZ':
        for u, v, data in out_edges:
            if data["pauli"] == pauliString[0]:
                assert 0 < data['weight'] <= 1.000001, f'data[weight] == {data["weight"]} > 1.0'
                ret = compute_covered_prob(graph, v, pauliString[1:])
                return data["weight"] * ret
        return d.Decimal(0)  # cannot be covered

    elif pauliString[0] == "I":
        _s = d.Decimal(0)
        for u, v, data in out_edges:
            assert 0 < data['weight'] <= 1.0, f'data[weight] == {data["weight"]} > 1.0'
            _s += data["weight"] * compute_covered_prob(graph, v, pauliString[1:])
        assert 0 <= _s <= 1.000001, f'0 <= {_s} <= 1.0 with node={curr_node} ({list(out_edges)}) and pauliString={pauliString}'
        return _s
    else:
        raise ValueError('Unknown Pauli operator "{}" in {}'.format(pauliString[0], pauliString))


def count_covering_paths(graph: nx.MultiDiGraph, curr_node: Hashable, pauli_str: PauliTerm) -> int:
    """count the number of paths covering pauliString str"""

    out_edges = graph.out_edges(curr_node, data=True)

    if out_edges:
        if pauli_str[0] in 'XYZ':
            for u, v, data in out_edges:
                if data["pauli"] == pauli_str[0]:
                    return count_covering_paths(graph, v, pauli_str[1:])  # find the rest of the probability
            return 0  # cannot be covered

        elif pauli_str[0] == "I":
            _s = 0
            for u, v, data in out_edges:
                _s += count_covering_paths(graph, v, pauli_str[1:])
            return _s
        else:
            raise ValueError('Unknown Pauli operator "{}"'.format(pauli_str[0]))

    else:  # no out edges
        if len(pauli_str) > 0 and ("X" in pauli_str or "Y" in pauli_str or "Z" in pauli_str):
            raise ValueError('No edges but pauli string not consumed')
        return 1


def get_covering_paths(graph: nx.MultiDiGraph, current_node: Hashable, pauli_str: PauliTerm, path: str, paths: Set[str]) -> Set[str]:
    out_edges = graph.out_edges(current_node, data=True)

    if out_edges:
        for u, v, data in out_edges:
            if pauli_str[0] == data["pauli"] or pauli_str[0] == 'I':
                new_path = path + data['pauli']
                paths = get_covering_paths(graph, v, pauli_str[1:], new_path, paths)
    else:
        paths.add(path)

    return paths


def does_path_go_through_edge(graph: nx.MultiDiGraph, pauli_str: PauliTerm, target_node: Hashable, target_pauli_op: PauliOp) -> bool:
    assert all(p in 'XYZ' for p in pauli_str)

    current_node = 0
    level = 0
    passed_edge = False
    while current_node != -1:
        out_edges = list(graph.out_edges(current_node, data=True))
        for u, v, data in out_edges:
            if data['pauli'] == pauli_str[level]:
                current_node = v
                level += 1
                if u == target_node and target_pauli_op == data['pauli']:
                    passed_edge = True
                break
        else:
            raise RuntimeError(f'Path {pauli_str} not found in graph.')

    return passed_edge


def get_covering_paths_through_edge(graph: nx.MultiDiGraph, pauli_str: PauliTerm, target_node: Hashable, target_pauli_op: PauliOp) -> Set[PauliTerm]:
    possible_paths = get_covering_paths(graph, 0, pauli_str, '', set())
    return set(p for p in possible_paths if does_path_go_through_edge(graph, p, target_node, target_pauli_op))


def count_covering_paths_through_edge(graph: nx.MultiDiGraph, pauli_str: PauliTerm, target_node: Hashable, target_pauli_op: PauliOp) -> int:
    possible_paths = get_covering_paths(graph, 0, pauli_str, '', set())
    return sum(1 for p in possible_paths if does_path_go_through_edge(graph, p, target_node, target_pauli_op))


def compute_covered_prob_through_edge(graph: nx.MultiDiGraph, pauli_str: PauliTerm, target_node: Hashable, target_pauli_op: PauliOp) -> d.Decimal:
    possible_paths = get_covering_paths(graph, 0, pauli_str, '', set())
    #print(f'{pauli_str} has {possible_paths}')
    filtered_paths = {p for p in possible_paths if does_path_go_through_edge(graph, p, target_node, target_pauli_op)}
    #print(f'{pauli_str} has {filtered_paths} through {target_node}_{target_pauli_op}')
    return sum(compute_covered_prob(graph, 0, p) for p in filtered_paths)

