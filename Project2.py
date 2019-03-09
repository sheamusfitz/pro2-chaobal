'''
Project 2 "Chaos Balls"
'''

import numpy as np
import matplotlib.pyplot as plyt

#	#	#	#	#	#	#	#	#	#	#
#		INITIAL CONDITIONS HERE			#
ma_raw = 1.0
mb_raw = 9.5

#				make xa0<xb0			#
xa0_raw = 1
xb0_raw = 3

va0_raw = 0
vb0_raw = 0

n_cols=6000 #the number of collisions to calculate

#Exact calculations: timesteps
dt = 0.1
part2steps = 4000

#higher values of tau give fewer points in the autocorrelation function
tau = 10


#										#
#	#	#	#	#	#	#	#	#	#	#



#	DON'T CHANGE DEFINITIONS AFTER HERE	#

if xa0_raw>xb0_raw:
	print('1-ERROR1 you did positions wrong, make xa>xb')

g_raw=9.81
e_raw=.5*ma_raw*va0_raw**2+.5*mb_raw*vb0_raw**2+(ma_raw*xa0_raw+mb_raw*xb0_raw)*g_raw
#print(f'2-total energy e_raw = {e_raw}')



def normalize(ma_raw,mb_raw,xa0_raw,xb0_raw,va0_raw,vb0_raw):
	g_raw=9.81
	e_raw=.5*ma_raw*va0_raw**2+.5*mb_raw*vb0_raw**2+(ma_raw*xa0_raw+mb_raw*xb0_raw)*g_raw
	#print(f'2-total energy e_raw = {e_raw}')
	# definition of normalization units:
	m=ma_raw+mb_raw
	x_n=e_raw/(m*g_raw)
	t_n=(e_raw/(m*g_raw**2))**.5
	v_n=(e_raw/m)**.5

	#normalizing raw input'd values
	ma=ma_raw/m
	mb=mb_raw/m
	
	xa0=xa0_raw/x_n
	xb0=xb0_raw/x_n
	
	va0=va0_raw/v_n
	vb0=vb0_raw/v_n
	#print(va0)
	
	#e_test = .5*ma*va0**2+.5*mb*vb0**2+ma*xa0+mb*xb0
	#print(f'total normalized energy \'{e_test}\' should be basically 1')
	
	return ma,mb,xa0,xb0,va0,vb0



t=np.zeros(n_cols)
xa=np.zeros(n_cols)
xb=np.zeros(n_cols)
va=np.zeros(n_cols)
vb=np.zeros(n_cols)


'''
-this is the list of ground bounce times. The point of this is so that I can
 track the full list, to be used for exact motion tracking. Similarly, there 
 are xigg and vigg for i==a,b
'''

tgg =np.zeros(part2steps)
xagg=np.zeros(part2steps)
xbgg=np.zeros(part2steps)
vagg=np.zeros(part2steps)
vbgg=np.zeros(part2steps)



def part1():
	global tgg,xagg,xbgg,vagg,vbgg
	#for some reason, THIS ONE needs to be declared global, but none of the others. go figure.
	
	t_c = 0
	
	#	normalizing initial conditions

	norm_vars=normalize(ma_raw,mb_raw,xa0_raw,xb0_raw,va0_raw,vb0_raw)
	ma=norm_vars[0]
	mb=norm_vars[1]
	
	t[0]=0
	
	xa[0]=norm_vars[2]
	xb[0]=norm_vars[3]
	
	tgg[0] = 0
	
	va[0]=norm_vars[4]
	vb[0]=norm_vars[5]
	
	'''
	t[i] is the time of the i'th collision
	t_t is short for 'time_temporary', and is always (I think) used for the time when _a hits the ground after the i'th collision
	xas and xbs are the bounce positions ie the locations of the two balls when _a hits the ground
	I hope you can guess what va_t and vb_t are
	'''
	
	k=1
	
	if va[0]==vb[0]:
		t_c = -1.0
		'''
		it isn't #actually# -1, but like... they won't collide at all if their initial velocities are the same. (they're moving in parallel in an accelerating reference frame, or however you want to phrase it)
		'''
	else:
		t_c = (t[0]*(va[0]-vb[0])-xa[0]+xb[0])/(va[0]-vb[0])
	
	for i in range(n_cols-1):
#		if i%modcount==0:
# 			print('i:',i)

		t_t = t[i]
		t_g = t[i]+va[i]+(va[i]**2+2*xa[i])**.5
		
		#print('5-',t[i],t_c)
		
		xa_t=xa[i]
		xb_t=xb[i]
		va_t=va[i]
		vb_t=vb[i]
		
		while t_g<t_c or t_c<=t_t:
			
			
			#print('entered the loop for i =',i)
			#now we have new temporary values
			xa_t = 0
			xb_t = xb_t+vb_t*(t_g-t_t)-.5*(t_g-t_t)**2
			va_t = (-1)*(va_t-(t_g-t_t))
			vb_t = vb_t-(t_g-t_t)
			
			if k<part2steps:
				tgg[k] = t_g
				xagg[k] = xa_t
				xbgg[k] = xb_t
				vagg[k] = va_t
				vbgg[k] = vb_t
				k+=1

			
			
			#our new "initial" time is our last t_g
			t_t = t_g
			
			#the next collision time is
			#t_c = (t_g*(va_t-vb_t)-xa_t+xb_t)/(va_t-vb_t)
			t_c = (t_g * (va_t - vb_t)-xa_t+xb_t)/(va_t-vb_t)
						
			#the next ground bounce time is
			t_g = t_g+va_t+(va_t**2+2*xa_t)**.5			
			
		if t_g==t_c:
			print('6-ERROR2, wow, you did the thing, weird')
			return
		
		if t_c<t_g and t_c>t_t:
			#print('yes')
			1
		else:
			print('NOOOO\n',t_t,t_c,t_g)
			return
		
		#now evolve one step
		t[i+1] = t_c
		
		xa[i+1] = xa_t + va_t*(t[i+1]-t_t) - .5*(t[i+1]-t_t)**2
		#print(xa[i+1])
		xb[i+1] = xb_t + vb_t*(t[i+1]-t_t) - .5*(t[i+1]-t_t)**2
#		print(xa[i+1],xb[i+1])
		xb[i+1] = xa[i+1]
		
		va_before = va_t - (t[i+1]-t_t)
		vb_before = vb_t - (t[i+1]-t_t)
		va[i+1]=((ma-mb)*va_before+2*mb*vb_before)/(ma+mb)
		vb[i+1]=((mb-ma)*vb_before+2*ma*va_before)/(ma+mb)
		
		#print('10-',t[i],t_c)
	
# 	tdif = np.zeros(n_cols-2)
# 	for i in range(n_cols-2):
# 		tdif[i] = t[i+1]-t[i]
# 		if tdif[i]>3:
# 			print(tdif[i])
#	print(np.amin(tdif))
# 	print(np.sqrt(np.mean(np.square(tdif))))
	
	#poincaré section==> psec

	'''
	plyt.figure()	
	plyt.plot(xb,vb,'k.', markersize=1)
	plyt.title('yeh')
	plyt.show()
	'''


#print(tgg.size,n_cols)

#print(xa)

t_ex  = np.zeros(part2steps)
xa_ex = np.zeros(part2steps)
xb_ex = np.zeros(part2steps)



def part2():
	#looping indeces
	i = 0
	j = 0
	n = 1
	
	t_ref = 0
	
	xa_ref=xa[0]
	xb_ref=xb[0]
	va_ref=va[0]
	vb_ref=vb[0]
	
	i_max = t.size-1
	j_max = tgg.size-1
	
	
	
	'''
	commenting b/c i need these in the next part as well, so they're getting initialized
	outside of part2()	
	'''
#	t_ex  = np.zeros(part2steps)
# 	xa_ex = np.zeros(part2steps)
# 	xb_ex = np.zeros(part2steps)
	xa_ex[0] = xa[0]
	xb_ex[0] = xb[0]
	va_ex = np.zeros(part2steps)
	vb_ex = np.zeros(part2steps)
	
	
	
	while n<part2steps and i<i_max and j<j_max:
		#print(t[i+1],tgg[j+1])
		
		if t[i+1] <= tgg[j+1]:
			#print('yeah')
			while n*dt < t[i+1] and n<part2steps:
				delta_t = n*dt-t_ref
				xa_ex[n] = xa_ref + va_ref*delta_t - 0.5*delta_t**2
				xb_ex[n] = xb_ref + vb_ref*delta_t - 0.5*delta_t**2
				n+=1
			i+=1
#
#			print('iii',i,n)
			t_ref = t[i]
			xa_ref = xa[i]
			xb_ref = xb[i]
			va_ref = va[i]
			vb_ref = vb[i]
		else:
			#print('else')
			while n*dt < tgg[j+1] and n<part2steps:
				delta_t = n*dt-t_ref
				xa_ex[n] = xa_ref + va_ref*delta_t - 0.5*delta_t**2
				xb_ex[n] = xb_ref + vb_ref*delta_t - 0.5*delta_t**2
				n+=1
			j+=1
#			print('jjjjjjjjj',j,n)

			t_ref = tgg[j]
			xa_ref = xagg[j]
			xb_ref = xbgg[j]
			va_ref = vagg[j]
			vb_ref = vbgg[j]
	
#	print(i_max,j_max)
	
	n=0
	for n in range(part2steps):
		t_ex[n]=n*dt
#		print(t_ex[n])
	'''	
	plyt.figure()
	plyt.plot(t_ex,xa_ex,t_ex,xb_ex,markersize=1)	
#	plyt.plot(t_ex,xa_ex,'.',t_ex,xb_ex,'.',markersize=1)
	#Axes.set_xlim(right=100)
	plyt.title('yeh')
	plyt.show()
#	print(xb_ex)
	'''
	


j_max = int(np.floor(part2steps/tau))-2
t_cor = np.zeros(j_max)
xa_cor= np.zeros(j_max)
xb_cor= np.zeros(j_max)


def autocor():
#	num_tau_steps is free to be changed.
# 	tau = 7
# 	
# 	j_max = int(np.floor(part2steps/tau))-2
# 	
# 	t_cor = np.zeros(j_max)
# 	xa_cor= np.zeros(j_max)
# 	xb_cor= np.zeros(j_max)
	
	global xa_cor,xb_cor
	
	j=0
	i=0
	
	for j in range(j_max):
		t_cor[j] = j*tau*dt
		i_max = (part2steps)-tau*(j+1)
		for i in range (i_max):
			xa_cor[j] += xa_ex[i]*xa_ex[i+tau*(j+1)]/i_max
			xb_cor[j] += xb_ex[i]*xb_ex[i+tau*(j+1)]/i_max
	
	'''
	-okay this seems to work, i have to normalize this ish now
	-i'm sure numpy has more efficient ways of doing this but i don't really 
	need to have a hyper efficient method here. 
	'''
	xa_cor_avg = np.mean(xa_cor)
	xb_cor_avg = np.mean(xb_cor)
	# for j in range(j_max):
# 		xa_cor[j] = xa_cor[j]-xa_cor_avg
# 		xb_cor[j] = xb_cor[j]-xb_cor_avg
	xa_cor = xa_cor-xa_cor_avg
	xb_cor = xb_cor-xb_cor_avg	
	xa_cor = xa_cor/np.sqrt(np.mean(np.square(xa_cor)))
	xb_cor = xb_cor/np.sqrt(np.mean(np.square(xb_cor)))

	
	
	
	
#	print(xb_ex)

'''		
part1()
part2()
autocor()
'''



#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	
#question 5
'''
ma_raw = 1.0
mb_raw = 9.7
xa0_raw = 1
xb0_raw = 3

va0_raw = 0
vb0_raw = 0

n_cols=10000 #the number of collisions to calculate
modcount = int(n_cols/20)

#Exact calculations: timesteps
dt = 0.1
part2steps = 4000

#higher values of tau give fewer points in the autocorrelation function
tau = 10

part1()
part2()
autocor()

plyt.figure()
plyt.subplot(311)
plyt.title(mb_raw)
plyt.plot(xb,vb,'k.', markersize=1)
plyt.subplot(312)
#	plyt.plot(t_ex,xa_ex,t_ex,xb_ex,markersize=1)	
plyt.plot(t_ex,xa_ex,'b.',t_ex,xb_ex,'r.',markersize=1)
#Axes.set_xlim(right=100)
#	print(xb_ex)
plyt.subplot(313)
plyt.plot(t_cor,xa_cor,'b.',t_cor,xb_cor,'r.',markersize=1)
#	plyt.plot(t_cor,xa_cor,t_cor,xb_cor)
plyt.show()
'''


#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	
# question 6

ma_raw = 1.0
mb_raw = 9
xa0_raw = 1
xb0_raw = 3

va0_raw = 0
vb0_raw = 0

#Exact calculations: timesteps
dt = 0.1
part2steps = 4000

#higher values of tau give fewer points in the autocorrelation function
tau = 5



plyt.figure()

some_step = 1.2
n = 12
for i in range(n):
	some_color = i/n
	#print(i)
	part1()
	ired   = 1-i/n
	igreen = 0
	iblue  = (1-ired) #.4*np.cos(3.14*(i+4)/12)+.5
	plyt.plot(xb,vb,'.', markersize=1,c=[ired, igreen, iblue])
	va0_raw += some_step
	#print(va0_raw)

plyt.xlim(right = 1)
plyt.show()



	