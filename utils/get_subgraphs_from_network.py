import sys
sys.setrecursionlimit(5000)

def get_cliques(size):
    cliques = set()
    if size == 3:
        for node in nodes:
            neigh = adjList[node]
            for i in range(len(neigh)):
                for j in range(i+1, len(neigh)):
                    u = neigh[i]
                    v = neigh[j]
                    if v in adjList[u]:
                        clique = [node, u, v]
                        clique.sort()
                        cliques.add(tuple(clique))

    elif size == 4:
        for node in nodes:
            neigh = adjList[node]
            for i in range(len(neigh)):
                for j in range(i+1, len(neigh)):
                    u = neigh[i]
                    v = neigh[j]
                    if not v in adjList[u]:
                        continue
                    for k in range(j+1, len(neigh)):
                        w = neigh[k]
                        if (w in adjList[v]) and (w in adjList[u]):
                            clique = [node, u, v, w]
                            clique.sort()
                            cliques.add(tuple(clique))

    elif size == 5:
        for node in nodes:
            neigh = adjList[node]
            for i in range(len(neigh)):
                for j in range(i+1, len(neigh)):
                    u = neigh[i]
                    v = neigh[j]
                    if not v in adjList[u]:
                        continue
                    for k in range(j+1, len(neigh)):
                        w = neigh[k]
                        if not ((w in adjList[v]) and (w in adjList[u])):
                            continue;
                        for l in range(k+1, len(neigh)):
                            z = neigh[l]
                            if (z in adjList[u]) and (z in adjList[v]) and (z in adjList[w]):
                                clique = [node, u, v, w, z]
                                clique.sort()
                                cliques.add(tuple(clique))

    else:
        print("Not supported")
    return cliques;


def get_subgraphs(size):
    subgrs = set()
    if size == 3:
        for node in nodes:
            neigh = adjList[node]
            for i in range(len(neigh)):
                for j in range(i+1, len(neigh)):
                    u = neigh[i]
                    v = neigh[j]
                    subgr = [node, u, v]
                    subgr.sort()
                    subgrs.add(tuple(subgr))

    elif size == 4:
        for node in nodes:
            neigh = adjList[node]
            for i in range(len(neigh)):
                for j in range(i+1, len(neigh)):
                    u = neigh[i]
                    v = neigh[j]
                    for k in range(j+1, len(neigh)):
                        w = neigh[k]
                        subgr = [node, u, v, w]
                        subgr.sort()
                        subgrs.add(tuple(subgr))

    else:
        return esu([], nodes, [], size)
    return subgrs


def esu(cur, left, neigh, remDepth):
    if remDepth == 0:
        return [cur]

    cand = []
    if len(cur) == 0:
        cand = left
    else:
        cand = list(set(neigh) & set(left))

    #print(cur, remDepth, neigh)
    if len(cand) != 0:
        v = cand[0]
        left2 = [x for x in left if x != v]
        out = []
        out = out + esu(cur, left2, neigh, remDepth)
        out = out + esu(cur+[v], left2, neigh+adjList[v], remDepth-1)
        return out

    return []







filen = sys.argv[1]
outn = sys.argv[2]
sub_size = int(sys.argv[3])
only_cliques = int(sys.argv[4])

tmpList = dict()
adjList = dict()
nodes = []

with open(filen) as file:
    for line in file:
        line = line.rstrip()
        token = line.split()
        token[0] = token[0].rstrip()
        token[1] = token[1].rstrip()
        if token[0] == token[1]:
            continue
        if token[0] not in tmpList:
            tmpList[token[0]] = set()
        if token[1] not in tmpList:
            tmpList[token[1]] = set()
        tmpList[token[0]].add(token[1])
        tmpList[token[1]].add(token[0])
    for key in tmpList:
        nodes.append(key)
        lst = list(tmpList[key])
        adjList[key] = lst
        #print(key, lst)

print(len(nodes))




if only_cliques:
    subgraphs = get_cliques(sub_size)
else:
    subgraphs = get_subgraphs(sub_size)

print(len(subgraphs))
with open(outn, 'w') as file:
    for subgr in subgraphs:
        for node in subgr:
            file.write(str(node)+" ")
        file.write("\n")
