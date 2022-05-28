import itertools
import os, shlex
import numpy as np
import numpy.typing as npt
from scipy.spatial import HalfspaceIntersection
from scipy.optimize import linprog

def nfg_parse(file_path):

    def list_parse(parts):
        assert(parts[0] == '{')
        results, parts = [], parts[1:]

        head, parts = parts[0], parts[1:]
        while head != '}':
            results.append(head)
            head, parts = parts[0], parts[1:]
        return parts, results

    with open(file_path, 'r') as file:
        parts = shlex.split(file.read())

        assert(parts[0] == 'NFG')
        assert(parts[1] == '1')
        assert(parts[2] == 'R')
        # parts[3] are comments
        assert(parts[4] == '{')

        parts, players    = list_parse(parts[4:])
        parts, strategies = list_parse(parts)

        strategies = list(map(int, strategies))
        print(strategies)
        outcomes = list(itertools.product(*[range(i) for i in reversed(strategies)]))
        outcomes = list(map(lambda x:tuple(reversed(x)), outcomes))
        print(outcomes)
        payoff = {}
        for outcome in outcomes:
            for player in range(len(strategies)):
                payoff[outcome + (player, )], parts = int(parts[0]), parts[1:]
        
        assert(parts == [])

        return strategies, payoff

def ne_output(file_path, results):
    with open(file_path, 'w') as file:
        for r in results:
            file.write(f"{r}\n")

def pure_nash(strategies, payoff):

    def get_payoff(outcome, player):
        return payoff[tuple(outcome) + (player, )]

    players = len(strategies)
    assert(players >= 3)

    strategies = [list(range(i)) for i in strategies]   # expansion of strategies list
    outcomes   =  list(itertools.product(*strategies))  # cartesian product of a series of lists
    results    = []

    # print(f"players: {players}")

    for outcome in outcomes:
        flag = True
        for player in range(players):
            for other_response in strategies[player]:
                other_outcome = list(outcome)
                other_outcome[player] = other_response

                if get_payoff(other_outcome, player) > get_payoff(outcome, player):
                    flag = False
                    break
            
            if flag == False:
                break
        
        if flag == True:
            result = []
            for player in range(players):
                tmp = [0] * len(strategies[player])
                tmp[outcome[player]] = 1
                result += tmp
            results.append(result)

    return list(map(list, results))

def mixed_nash(strategies, payoff):

    def powerset(iterable):
        "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        s = list(iterable)
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))
    
    def build_polytope(M: npt.NDArray) -> npt.NDArray:
        n, m = M.shape
        b = np.append(-np.ones(n), np.zeros(m))
        M = np.append(M, -np.eye(m), axis = 0)
        return np.column_stack((M, b.transpose()))
    
    def find_points(M: npt.NDArray):
        nvec = np.reshape(
            np.linalg.norm(M[:, :-1], axis = 1),
            (M.shape[0], 1)
        )
        c = np.zeros((M.shape[1], ))
        c[-1] = -1
        A = np.hstack((M[:, :-1], nvec))
        b = -M[:, -1:]
        res = linprog(c, A_ub = A, b_ub = b)
        return res.x[:-1]
    
    def get_label(v: npt.NDArray, tM: npt.NDArray):
        b = tM[:, -1]
        M = tM[:, :-1]
        return set(np.where(np.isclose(np.dot(M, v), -b))[0])
    
    def vertices(M: npt.NDArray):
        points = find_points(M)
        hs = HalfspaceIntersection(M, points)
        hs.close()
        results = []
        for v in hs.intersections:
            if not np.all(np.isclose(v, 0)) and max(v) < np.inf:
                results.append((v, get_label(v, M)))
        return results

    players = len(strategies)
    assert(players == 2)

    outcomes = list(itertools.product(*[range(i) for i in strategies]))

    n, m = strategies[0], strategies[1]

    A = [[payoff[(a, b, 0)] for b in range(m)]
            for a in range(n)]
    B = [[payoff[(a, b, 1)] for b in range(m)]
            for a in range(n)]
    
    A = np.array(A)
    B = np.array(B)

    print(A)
    print(B)

    if np.min(A) < 0: A = A + abs(np.min(A))
    if np.min(B) < 0: B = B + abs(np.min(B))

    n, m = A.shape
    total = n + m
    labels = set(range(total))

    A_polytope = build_polytope(A)
    B_polytope = build_polytope(B.transpose())
    result = []

    for Bx, By in vertices(B_polytope):
        tmp = set((i + n) % total for i in By)

        for Ax, Ay in vertices(A_polytope):
            if tmp.union(Ay) == labels:
                arr = np.append(Bx / sum(Bx), Ax / sum(Ax))
                result.append(arr.tolist())

    return result

def nash(in_path, out_path):
    # load file
    strategies, payoff = nfg_parse(in_path)
    players = len(strategies)

    # get NE
    if players == 2:
        results = mixed_nash(strategies, payoff)
    else:
        results = pure_nash (strategies, payoff)
    
    # write file
    ne_output(out_path, results)

if __name__ == '__main__':
    for f in os.listdir('examples'):
        if f.endswith('.nfg'):
            print(f)
            nash('examples/'+f, 'output/'+f.replace('nfg','ne'))
    
    for f in os.listdir('examples'):
        if f.endswith('.ne'):
            print(f)
            