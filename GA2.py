import numpy as np
import random
import math

class DNA(object):
    def __init__(self, a, b, c, cross_rate, mutate_rate):
        self.a = a
        self.b = b
        self.c = c
        self.cross_rate = cross_rate
        self.mutate_rate = mutate_rate
        self.dna = []
        self.cons_num = len(self.b)     # number of constrains
        self.x_num = len(self.c)        # length of DNA
        for i in range(self.x_num):     # randomly initialize DNA
            self.dna.append(random.randint(0, 1))

    def get_value(self):
        value = 0
        for i in range(self.x_num):
            value = value + self.c[i] * self.dna[i]
        for i in range(self.cons_num):
            sum = 0     # numbers of constrains be broken
            for j in range(self.x_num):
                sum = sum + self.a[i][j] * self.dna[j]
            if sum > self.b[i]:
                value = value - 3.33    # when any constrain is broken, get negative reward
        return value

    def mutate(self):
        for i in range(self.x_num):
            if random.random() < self.mutate_rate:
                if self.dna[i] == 0:        # randomly mutate
                    self.dna[i] = 1
                else:
                    self.dna[i] = 0

    def cross_over(self, neigh):
        new_dna = []
        for i in range(self.x_num):
            new_dna.append(random.randint(0, 1))    # find the other parent
        for i in range(self.x_num):
            if random.random() < 0.5:       #  50% probability inherit from each parent
                new_dna[i] = self.dna[i]
            else:
                new_dna[i] = neigh.dna[i]
        return new_dna

class Population(object):
    def __init__(self, a, b, c, cross_rate, mutate_rate, pop_size):
        self.a = a
        self.b = b
        self.c = c
        self.cross_rate = cross_rate        # crossover probability
        self.mutate_rate = mutate_rate      # mutate probability
        self.cons_num = len(self.b)
        self.x_num = len(self.c)
        self.pop_size = pop_size
        self.DNA_list = []          # population
        for i in range(self.pop_size):
            temp_DNA = DNA(a, b, c, self.cross_rate, self.mutate_rate)
            self.DNA_list.append(temp_DNA)
        self.fitness = []

    def get_value(self):
        self.fitness = []       # update fitness of all DNAs
        for i in range(self.pop_size):
            self.fitness.append(self.DNA_list[i].get_value())

    def cross_over(self):
        for i in range(self.pop_size):
            if random.random() < self.cross_rate:   # happens
                neigh = random.randint(0, self.pop_size - 1)
                while(neigh == i):
                    neigh = random.randint(0, self.pop_size - 1)    # choose another parent
                self.DNA_list[i].dna = self.DNA_list[i].cross_over(self.DNA_list[neigh])

    def mutate(self):
        for i in range(self.pop_size):
            self.DNA_list[i].mutate()

    def select(self):
        new_dna_list = []
        temp_fit = self.fitness-np.min(self.fitness)+0.01
        index = np.random.choice(a=self.pop_size, size=self.pop_size, replace=True, p=temp_fit/np.sum(temp_fit))
        for i in range(len(index)):
            new_dna_list.append(self.DNA_list[index[i]])
        self.DNA_list = new_dna_list

    def get_opt_dna(self):
        arr_aa = np.array(self.fitness)
        index = np.argmax(arr_aa)
        return np.max(arr_aa), self.DNA_list[index].dna

    def dec_to_bin(self, sb):
        my_list = []
        for i in range(len(c)):
            my_list.append(0)
        temp = sb
        for i in range(len(c)):
            check = math.pow(2, len(c) - i - 1)
            if temp < check:
                my_list[i] = 0
            else:
                my_list[i] = 1
                temp = temp - check
        return my_list

    def ana_sol(self):
        b_fitness = -1000000
        b_dna = []
        for i in range(int(math.pow(2, self.DNA_list[0].x_num))):
            my_list = self.dec_to_bin(i)
            value = 0
            for i in range(self.x_num):
                value = value + self.c[i] * my_list[i]
            for i in range(self.cons_num):
                sum = 0
                for j in range(self.x_num):
                    sum = sum + self.a[i][j] * my_list[j]
                if sum > self.b[i]:
                    value = value - 3.33
            if value>b_fitness:
                b_fitness = value
                b_dna = my_list
        return b_fitness, b_dna

    def print_DNA(self):
        for i in range(self.pop_size):
            print('DNA', i, self.DNA_list[i].dna)

    def inject_best(self, g_best_fit, g_best_dna):
        arr_aa = np.array(self.fitness)
        index = np.argmin(arr_aa)
        self.DNA_list[index].dna = list(g_best_dna)
        self.fitness[index] = g_best_fit

# GA
opt_iter = 50        # iteration
pop_size = 10        # number of DNA
cross_rate = 0.3        # crossover probability
mutate_rate = 0.01      # mutate probability

a = [[2, -6, 3, -4, -1, 2], [5, 3, -1, -3, 2, -1], [-5, 1, -4, 2, -2, 1]]
b = [-2, 2, -3]
c = [3, 5, 6, 9, 10, 10]

pop = Population(a, b, c, cross_rate, mutate_rate, pop_size)
g_best_fit = -10000
g_best_dna = []
for i in range(opt_iter):
    pop.cross_over()
    pop.mutate()
    pop.get_value()
    best_fit, best_dna = pop.get_opt_dna()
    if best_fit > g_best_fit:       # record global best
        # print('update')
        g_best_fit = best_fit       # record global best
        g_best_dna = best_dna
    print('iter', i, 'fitness', g_best_fit, 'solution', g_best_dna)
    # pop.print_DNA()
    pop.inject_best(g_best_fit, g_best_dna)
    pop.select()
    pop.get_value()
b_fitness, b_dna = pop.ana_sol()
print('analytical fitness', b_fitness, 'solution', b_dna)