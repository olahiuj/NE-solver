import itertools
import os
import shlex

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
        outcomes = list(itertools.product(*[range(i) for i in strategies]))
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
    players = len(strategies)
    assert(players == 2)
    return []

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

    for f in os.listdir('input'):
        if f.endswith('.nfg'):
            print(f)
            nash('input/'+f, 'output/'+f.replace('nfg','ne'))