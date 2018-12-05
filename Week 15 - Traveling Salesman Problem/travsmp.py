
import numpy as np
import matplotlib.pyplot as plt

L = 15 
T0 = L/8.0 # intial temperature as stated in the book

# randomly scatter cities (True) in an LxL lattice
canvass = np.random.random([L,L])<0.5 

#number of cities
N = np.sum(canvass) 

#location of cities
cities = np.array(zip(np.where(canvass)[0],np.where(canvass)[1]))
distances = []

def distance(cities):
    '''Calculate the total distance traveled.'''
    
    d = 0.0
    city1,city2 = cities,np.roll(cities,-1,axis=0)
    x1,y1 = city1[:,1],city1[:,0]
    x2,y2 = city2[:,1],city2[:,0]
    d =  np.sqrt((x2-x1)**2+(y2-y1)**2) 
    return np.sum(d)

def otherroutes(cities):
    '''Check out other routes'''
    
    new_cities = np.copy(cities)
    p = np.random.randint(N)
    l = np.random.randint(N//2)
    new_cities[p:p+l+1] = cities[p:p+l+1][::-1]
    return new_cities

def iterate(cities,T):
    new_cities = otherroutes(cities)    
    if distance(new_cities)<distance(cities):
        cities = new_cities    
    return cities
    
temp = []
T = T0
for i in xrange(1500):
    temp.append(T)
    cities = iterate(cities,T)
    distances.append(distance(cities))
    if i%10==0: 
        T = T*0.95
        
plt.plot(temp)  
plt.ylabel('T',size='large')
plt.xlabel('Iterations',size='large')
plt.show()
plt.close()

plt.plot(distances)
plt.ylabel('Distance traveled',size='large')
plt.xlabel('Iterations',size='large')
plt.show()
plt.close()

#route of the traveling salesman
x = np.array(list(cities[:,1]) + list(cities[-1]))
y = np.array(list(cities[:,0])  + list(cities[-1]))
plt.scatter(x,y,clip_on=0,marker='.')
plt.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], scale_units='xy', angles='xy', scale=1,width=0.003,headwidth=10,color='red')
plt.title('distance traveled = %.3f \n %i cities'%(distance(cities),N))
plt.axis('off')
plt.show()