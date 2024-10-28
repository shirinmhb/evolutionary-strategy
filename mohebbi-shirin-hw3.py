import numpy as np
import random
import matplotlib.pyplot as plt

class Es():
  def __init__(self, m, w, wv2, vAll, cAll, wAll, b, alpha, numPop, epsilon, obj):
    self.w = w
    self.wv2 = wv2
    self.m = m
    self.vAll = vAll
    self.cAll = cAll
    self.wAll = wAll
    self.alpha = alpha
    self.b = b
    self.numPop = numPop
    self.epsilon = epsilon
    self.obj = obj
    self.avgFitnesses = []
    self.maxFitness = []

  def obj1(self, r, n):
    R = []
    for i in range(len(r)):
      R.append(1 - (1 - r[i])**n[i])
    res = (R[0] * R[1]) + (R[2] * R[3]) + (R[0] * R[3] * R[4]) + (R[1] * R[2] * R[4]) - (R[0] * R[1] * R[2] * R[3]) - (R[0] * R[1] * R[2] * R[4]) - (R[0] * R[1] * R[3] * R[4]) - (R[0] * R[2] * R[3] * R[4]) - (R[1] * R[2] * R[3] * R[4]) + (2 * R[0] * R[1] * R[2] * R[3] * R[4])  
    return res
  

  def obj2(self, r, n):
    s = 1
    for i in range(self.m):
      a = (1 - (1 - r[i])**n[i])
      s *= a
    return s

  def constraint1(self, n):
    s = 0
    for i in range(self.m):
      s += (self.wv2[i] * n[i]**2)
    if s - self.vAll <= 0:
      return True
    return False

  def constraint2(self, r, n):
    s = 0
    for i in range(self.m):
      a = self.alpha[i] * (float(-1000)/np.log(r[i]))**self.b[i] * (n[i] + np.exp(0.25 * n[i]))
      s += a
    if s - self.cAll <= 0:
      return True
    return False

  def constraint3(self, n):
    s = 0
    for i in range(self.m):
      s += (self.w[i] * n[i] * np.exp(0.25 * n[i]))
    if s - self.wAll <= 0:
      return True
    return False

  def initialization(self, m):
    #each chrom is like this: [[r1,r2,r3,r3,r5],[sr1,sr2,sr3,sr4,sr5],[n1,n2,n3,n4,n5],[sn1,sn2,sn3,sn4,sn5]]
    pop = []
    for i in range(self.numPop):
      r = []
      sr = []
      n = []
      sn = []
      for j in range(5):
        r.append(random.uniform(0, 1))
        n.append(random.randint(1, m))
        sr.append(random.uniform(-0.01, 0.01))
        sn.append(random.uniform(-1,1))
      pop.append([r, sr, n, sn])
    self.population = pop

  def fitness(self, chrom):
    c1 = self.constraint1(chrom[2])
    c2 = self.constraint2(chrom[0], chrom[2])
    c3 = self.constraint3(chrom[2])
    if (c1 and c2 and c3):
      if (self.obj == 1):
        return self.obj1(chrom[0], chrom[2])  #complex system
      else:
        return self.obj2(chrom[0], chrom[2])  #series system
    else:  #if any constraint voilated, with would be negative number 
      return -10

  def parentSelection(self):   #uniform random
    selectedParent = []
    for i in range(7 * self.numPop):
      p1 = self.population[random.randint(0, numPop - 1)]
      p2 = self.population[random.randint(0, numPop - 1)]
      selectedParent.append([p1,p2])
    self.selectedParent = selectedParent

  def applyCrossOver(self):
    children = []
    for parents in self.selectedParent:
      child = self.crossOver(parents)
      children.append(child)
    self.offSprings = children

  def crossOver(self, parents):   #local intermediary
    #[[r1,r2,r3,r4,r5],[sr1,sr2,sr3,sr4,sr5],[n1,n2,n3,n4,n5],[sn1,sn2,sn3,sn4,sn5]]
    r = []
    sr = []
    n = []
    sn = []
    for i in range(5):
      newR = ( parents[0][0][i] + parents[1][0][i] )/ 2
      if newR > 1:  #r should be between 0 and 1
        newR = 1
      elif newR < 0:
        newR = 0
      r.append(newR)
      newN = round( ( parents[0][2][i] + parents[1][2][i] ) / 2 )  #n should be integer number
      # if newN > self.m:
      #   newN = self.m
      if newN < 1:
        newN = 1
      n.append(newN) 
      sr.append(( parents[0][1][i] + parents[1][1][i] )/2)
      sn.append(( parents[0][3][i] + parents[1][3][i] )/2)
    child = [r, sr, n, sn]
    return child

  def applyMutation(self):
    for child in self.offSprings:
      self.mutation(child)

  def mutation(self, child):
    #mutate sigmas:
    overAllLearning = 0.47 * np.random.normal(0.0, 1.0) 
    for i in range(5):
      coordinateLearning = 0.22 * np.random.normal(0.0, 1.0)
      child[1][i] = child[1][i] * np.exp(overAllLearning + coordinateLearning)
      if (child[1][i] < self.epsilon):
        child[1][i] = self.epsilon
      coordinateLearning = 0.22 * np.random.normal(0.0, 1.0)
      child[3][i] = child[3][i] * np.exp(overAllLearning + coordinateLearning)
      if (child[3][i] < self.epsilon):
        child[3][i] = self.epsilon
    #mutate x:
    for i in range(5):
      child[0][i] = child[0][i] + (child[1][i] * np.random.normal(0.0, 1.0))
      if child[0][i] < 0:
        child[0][i] = 0
      elif child[0][i] > 1:
        child[0][i] = 1
      child[2][i] = round(child[2][i] + (child[3][i] * np.random.normal(0.0, 1.0)))
      if child[2][i] < 1:
        child[2][i] = 1
      # elif child[2][i] > self.m:
      #   child[2][i] = self.m

  def survivalSelection(self):
    childWithFitness = []
    for child in self.offSprings:
      f = self.fitness(child)
      # print ("f", f)
      childWithFitness.append((child, f))
    childWithFitness.sort(key=lambda tup: tup[1], reverse=True)
    newpop = []
    avgFitness = 0
    for i in range(numPop):
      newpop.append(childWithFitness[i][0])
      avgFitness += childWithFitness[i][1]
    self.avgFitnesses.append(avgFitness/self.numPop)
    self.maxFitness.append(childWithFitness[0][1])
    print(childWithFitness[0][0])
    self.population = newpop

  def plot(self, title):
    if len(self.maxFitness) == 1:
      plt.scatter(x=0, y=self.maxFitness[0], color = "#00ace6", label = "max fitness")
      plt.scatter(x=0, y=self.avgFitnesses[0], color = "#9933ff", label = "avg fitness")
      plt.xlabel('iteration') 
      plt.ylabel('fitness')
      plt.title(title) 
      plt.legend()
      plt.show()
      return
    plt.plot(self.iterations, self.avgFitnesses, color = "#9933ff", label = "avg fitness")
    plt.plot(self.iterations, self.maxFitness, color = "#00ace6", label = "max fitness")
    plt.xlabel('iteration') 
    plt.ylabel('fitness')
    plt.title(title) 
    plt.legend()
    print("max fitness", max(self.maxFitness))
    plt.show()

  def main(self):
    self.initialization(5)
    for i in range(50):
      self.parentSelection()
      self.applyCrossOver()
      self.applyMutation()
      self.survivalSelection()
    self.iterations = list(range(0, i+1))
    self.plot("")
    # print (self.population[0], self.fitness(self.population[0]))
  
w = [7, 8, 8, 6, 9]
wv2 = [1, 2, 3, 4, 2]
m = 5
vAll = 110
cAll = 175
wAll = 200
b = [1.5, 1.5, 1.5, 1.5, 1.5]
alpha = [2.330 * 10**-5, 1.450 * 10**-5, 0.541 * 10**-5, 8.050 * 10**-5, 1.950 * 10**-5]
numPop = 100
epsilon = 0.005
es = Es(m, w, wv2, vAll, cAll, wAll, b, alpha, numPop, epsilon, 1)
es.main()


# w = [7, 8, 8, 6, 9]
# wv2 = [1, 2, 3, 4, 2]
# m = 5
# vAll = 110
# cAll = 175
# wAll = 200
# b = [1.5, 1.5, 1.5, 1.5, 1.5]
# alpha = [2.330 * 10**-5, 1.450 * 10**-5, 0.541 * 10**-5, 8.050 * 10**-5, 1.950 * 10**-5]
# numPop = 100
# epsilon = 0.005
# es = Es(m, w, wv2, vAll, cAll, wAll, b, alpha, numPop, epsilon, 2)
# es.main()


