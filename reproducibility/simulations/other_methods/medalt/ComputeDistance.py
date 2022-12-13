import copy


def matrixbuilder(node):
    matrix = []
    for node1 in node:
        temp = []
        for node2 in node:
            temp.append(dist(node[node1], node[node2]))
        matrix.append(temp)
    return node.keys(), matrix


def dist(node1, node2):
    d = 0
    for i in range(0, len(node1)):
        d = d+ disthelper(node1[i], node2[i])
    return d

def disthelper(node1, node2):
    if 0 in node1 or 0 in node2:
        return zerodisthelper(node1, node2)
    return distcalc(node1, node2)


def distcalc(node1, node2):
    assert len(node1) == len(node2)
    if len(node1) == 1:
        return abs(node1[0] - node2[0])
    else:
        d = 0
        newlist = copy.deepcopy(node1)
        for i in range(0, len(node2)):
            newlist[i] -= node2[i]
        while newlist:
            if newlist[0] == 0:
                newlist.pop(0)
            elif newlist[0] > 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] > 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] -= 1
                d += 1
            elif newlist[0] < 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] < 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] += 1
                d += 1
        return abs(d)


def zerodisthelper(node1, node2):
    n1 = copy.deepcopy(node1)
    n2 = copy.deepcopy(node2)
    dist = 0
    temp1 = []
    temp2 = []
    while n1:
        x1 = n1.pop()
        x2 = n2.pop()
        if x1 == 0:
            if x2 == 0:
                temp1.append(x1)
                temp2.append(x2)
            else:
                return 1000000
        else:
            temp1.append(x1)
            temp2.append(x2)
    return distcalc(temp1, temp2)
