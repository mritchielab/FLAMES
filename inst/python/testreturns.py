
def t1(ls):
    d = {}
    i = 1
    for key in ls:
        d[key] = i
        i += 1
    return i, d

def t2(d):
    print(d)
    print(type(d))


def t3(f):
    for i in range(10):
        print(i)
    print(f)