import math

# Experimental values
a = 1.09*1000
p = 0.33
e = math.sqrt(14.4)

# Defining Function
def V(r):
	return -e**2/r + a*(e**(-r/p))


# Defining derivative of function
def f(r):
	return e**2/r**2 - (a/p) * e**(-r/p)


# Implementing Newton Raphson Method

def newtonRaphson(r0, error, N):
	global r1
	print('\n\n*** NEWTON RAPHSON METHOD IMPLEMENTATION ***')
	step = 1
	flag = 1
	condition = True
	while condition:
		if f(r0) == 0.0:
			print('Divide by zero error!')
			break
		
		r1 = r0 - V(r0) / f(r0)
		print('Iteration-%d, r1 = %0.6f Å and V(r1) = %0.6f' % (step, r1, V(r1)))
		r0 = r1
		step = step + 1
		
		if step > N:
			flag = 0
			break
		
		condition = abs(V(r1)) > error
	
	if flag == 1:
		print('\nRequired root is: %0.8f' % r1)
	else:
		print('\nNot Convergent.')

# Input Section
r0 = input('Enter Guess: ')
error = input('Tolerable Error: ')
N = input('Maximum Step: ')

# Converting r0 and e to float
r0 = float(r0)
error = float(error)

# Converting N to integer
N = int(N)

# Note: You can combine above three section like this
# r0 = float(input('Enter Guess: '))
# e = float(input('Tolerable Error: '))
# N = int(input('Maximum Step: '))

# Starting Newton Raphson Method
newtonRaphson(r0, error, N)