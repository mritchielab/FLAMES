from test2 import add
import pysam

def python_test(vector):
  print(vector)
  s = 0
  for i in range(len(vector)):
    s += vector[i]
  return s

def python_test2(x, y):
  c = add(x,y)

  print c
