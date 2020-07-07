import matplotlib.pyplot as plt
from math import sqrt, pi
import csv


def three_body_problem(M, x, M_j, planet):

	delta_t = 0.001
	n_steps = int(11860) 

	M_s = 1.99 * 10**30
	M_e = M
	
	# Initial positions and velocities
	x_e = x
	y_e = 0.0
	x_j = 5.2
	y_j = 0.0
	v_ex = 0.0
	v_ey = sqrt(4 * pi**2 / x_e)
	v_jx = 0.0
	v_jy = sqrt(4 * pi**2 / x_j)

	f = open(planet + '.csv', 'w+')

	for i in range(1, n_steps+1):

		# Save positions to file
		f.write(str(x_e) + ' ' + str(y_e) + ' ' + str(x_j) + ' ' + str(y_j) + '\n')

		# Calculate distances
		r_e = sqrt(x_e**2 + y_e**2)
		r_j = sqrt(x_j**2 + y_j**2)
		r_ej = sqrt((x_e - x_j)**2 + (y_e - y_j)**2)

		# Compute new velocities
		v_ex_next = v_ex - (4 * pi**2 * x_e)/(r_e**3) * delta_t - (4 * pi**2 * (M_j / M_s) * (x_e - x_j))/(r_ej**3) * delta_t
		v_ey_next = v_ey - (4 * pi**2 * y_e)/(r_e**3) * delta_t - (4 * pi**2 * (M_j / M_s) * (y_e - y_j))/(r_ej**3) * delta_t
		v_jx_next = v_jx - (4 * pi**2 * x_j)/(r_j**3) * delta_t - (4 * pi**2 * (M / M_s) * (x_j - x_e))/(r_ej**3) * delta_t
		v_jy_next = v_jy - (4 * pi**2 * y_j)/(r_j**3) * delta_t - (4 * pi**2 * (M / M_s) * (y_j - y_e))/(r_ej**3) * delta_t

		# Compute new positions with EC method
		x_e_next = x_e + v_ex_next * delta_t
		y_e_next = y_e + v_ey_next * delta_t
		x_j_next = x_j + v_jx_next * delta_t
		y_j_next = y_j + v_jy_next * delta_t

		x_e = x_e_next
		y_e = y_e_next
		x_j = x_j_next
		y_j = y_j_next
	
		v_ex = v_ex_next
		v_ey = v_ey_next
		v_jx = v_jx_next
		v_jy = v_jy_next

		print i, 'of', n_steps, 'steps done for', planet
	
	f.close()	


def plot_orbits(planet):

	f = open(planet + '.csv')
	reader = csv.reader(f, delimiter = ' ')
	data = [(float(c1), float(c2), float(c3), float(c4)) for c1, c2, c3, c4 in reader]

	# Both orbits
	plt.figure(1)
	plt.xlabel('x (AU)')
	plt.ylabel('y (AU)')
	plt.title(planet);	
	plt.grid(True)
	plt.axis('equal')
	for row in data:
		plt.plot(row[0], row[1], '.', color = 'b', markersize = 1)
		plt.plot(row[2], row[3], '.', color = 'r', markersize = 1)
	plt.savefig(planet + '_1.png')
	plt.clf()
	
	# Just smaller planet
	plt.figure(2)
	plt.xlabel('x (AU)')
	plt.ylabel('y (AU)')
	plt.title(planet);	
	plt.grid(True)
	plt.axis('equal')
	plt.xlim((-1.7, 1.7))
	plt.ylim((-1.7, 1.7))
	for row in data:
		plt.plot(row[0], row[1], '.', color = 'b', markersize = 1)
	plt.savefig(planet + '_2.png')		
	plt.clf()
	
	f.close()
		

def calc_eccentricitiy(planet):
	f = open(planet + '.csv')
	reader = csv.reader(f, delimiter = ' ')
	data = [(float(c1), float(c2), float(c3), float(c4)) for c1, c2, c3, c4 in reader]	

	# Earth/Mars
	x = [row[0] for row in data]
	y = [row[1] for row in data]
	a = (abs(min(x)) + max(x))/2
	b = (abs(min(y)) + max(y))/2
	if b < a:
		ecc = sqrt(1 - (b**2 / a**2))
	else:
		ecc = sqrt(1 - (a**2 / b**2))
	print planet, a, b, ecc
	
	# Jupiter
	x = [row[2] for row in data]
	y = [row[3] for row in data]
	a = (abs(min(x)) + max(x))/2
	b = (abs(min(y)) + max(y))/2
	if b < a:
		ecc = sqrt(1 - (b**2 / a**2))
	else:
		ecc = sqrt(1 - (a**2 / b**2))
	print 'Jupiter', a, b, ecc


planet = 'Earth'
three_body_problem(6.0 * 10**24, 1.0,  1.9 * 10**27, planet)
#plot_orbits(planet)
calc_eccentricitiy(planet)


#planet = 'Mars'
#three_body_problem(6.6 * 10**23, 1.52, 1.9 * 10**27,  planet)
#plot_orbits(planet)
#calc_eccentricitiy(planet)
	

	
	


