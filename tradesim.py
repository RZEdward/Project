import matplotlib.pyplot as plt
import numpy
import random

bal = 10000
R = 2
risk = 0.02
winrate = 0.5
sims = 10000
steps = 500 #approx num of trades taken in 5 years - 2 trades taken Mon - Fri

plt.figure(1)
count = 0
for j in range(sims):
    trade = [0]
    balance = [bal]
    for i in range(steps):
        trade.append(i)
        acc = balance[i] 
        if random.uniform(0,1) < winrate:
            acc += R*risk*acc
        else:
            acc -= risk*acc
        balance.append(acc)
    count += balance[steps] - bal
    plt.plot(trade,balance)
av = count/sims
result = (av*100)/bal
label = "Average Return = " + str(result) + "%"
plt.title(label)
plt.show()



    
