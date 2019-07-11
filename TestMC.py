#!/usr/bin/env python3
'''
	TestMC.py

	Written by Marcos Flores,  04-Jul-2019
'''

import numpy as np
import random as rdm

import matplotlib
import matplotlib.pyplot as plt


def randprams(minE,maxE):
	# Outputs a random energy (in MeV) and a random angle (b/w -pi to pi)
	return [rdm.uniform(minE,maxE),rdm.uniform(-1,1)]

def main():

	mpi = 134.9766 # MeV

	num = int(1e3)

	cangle = np.zeros(2*num)
	energies = np.zeros(2*num)
	pizeromom = np.zeros(2*num)

	E = np.linspace(mpi,mpi*1e4,num)

	rnum = np.zeros(num)


	for i in range(0,num):
	
		[Epi,ctheta] = randprams(mpi,mpi*1e4)

		gamma = Epi/mpi
		beta = np.sqrt(1-1/(gamma**2))

		Egamma1 = gamma*mpi/2*(1 + beta*ctheta)
		Egamma2 = gamma*mpi/2*(1 - beta*ctheta)

		pizeromom[2*i] = Epi**2 - mpi**2
		pizeromom[2*i+1] = Epi**2 - mpi**2
		cangle[2*i] = ctheta
		cangle[2*i+1] = ctheta
		
		rnum[i] = ctheta

		energies[2*i] = Egamma1
		energies[2*i + 1] = Egamma2

	# plt.plot(range(0,num),rnum)

	fig1 = plt.figure(1)

	plt.scatter(cangle,energies,0.2)
	plt.xlim([-1.5, 1.5])
	plt.xlabel(r"cos $\theta$")
	plt.ylabel(r"$\gamma$ Energy (MeV)")

	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	# fig1.show()

	fig2 = plt.figure(2)

	plt.hist(energies, bins='auto')
	plt.xlabel(r"$\gamma$ Energy (MeV)")
	plt.ylabel("Count (#)")
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

	plt.plot(E,2/E,color="red")

	fig2.show()

	fig3 = plt.figure(3)	

	plt.hist(cangle, bins='auto',density=True)
	plt.xlabel(r"$\cos\theta$")

	plt.ylabel("Count (#)")
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

	fig3.show()

	input()

if __name__ == "__main__":

	main()