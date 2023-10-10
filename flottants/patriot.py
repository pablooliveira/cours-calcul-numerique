import matplotlib.pyplot as plt

def patriot(n):
    s = 0.0
    for i in range(n):
        s += 0.1
    return s

def exact(n):
    return n/10.

def relative_error(n): return abs((patriot(n) - exact(n))/exact(n))


ns = [2**i-1 for i in range(1, 20)]
data = {
    'n': ns,
    'abs-error':  list(map(relative_error, ns))
}
plt.plot('n', 'abs-error', data=data)
plt.xlabel('n')
plt.ylabel('abs-error')
plt.show()
