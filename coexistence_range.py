

#!/usr/bin/python3

import tkinter

import math

def destruction():

	window.destroy()

	return()

def checkData():

	global stepD

	global stepy

	def errorDestruction():

		errorWindow.destroy()

		return()

	g1 = float(entryg1.get())

	g2 = float(entryg2.get())

	k1 = float(entryk1.get())

	k2 = float(entryk2.get())

	if (g1/g2)>=(k1/k2):

		errorWindow = tkinter.Tk()

		errorLabel = tkinter.Label(errorWindow, text = 'Error: ratio of growth rates g\u2081/g\u2082 is equal to or greater than ratio of Monod constants K\u2081/K\u2082.\nIt is impossible for species with these characteristics to coexist on the same food source.')

		errorLabel.grid(row = 0, column = 0)

		errorButton = tkinter.Button(errorWindow, text = 'OK', command = errorDestruction)

		errorButton.grid(row = 1, column = 0)

	elif (g2/g1)>=1:

		errorWindow = tkinter.Tk()

		errorLabel = tkinter.Label(errorWindow, text = 'Error: ratio of growth rates g\u2082/g\u2081 is equal to or greater than 1.\nConsider inverting order of species.')

		errorLabel.grid(row = 0, column = 0)

		errorButton = tkinter.Button(errorWindow, text = 'OK', command = errorDestruction)

		errorButton.grid(row = 1, column = 0)

	else:

		possible = []

		possibleD = []

		possibler0 = []

		gamma = g2/g1 

		k = k1/k2 

		y0 = (gamma*k-1)/((1-gamma)*k)

		y = (1+stepy)*y0

		ymax = 10*k

		while y<ymax:

			m = gamma*(k-1)/((1-gamma)*(y*k+1)-gamma*(k-1))

			Dmax = (y*k+1)**m

			D = 1+stepD*(Dmax-1)

			while D<=Dmax and D<11000:

				logD = math.log(D)

				less = (y*k*D/(D-1)+1)*logD/((y*k*D/(D-1)+k)*logD+(k-1)*math.log(y*k+1))

				more = ((y*k*D/(D-1)+1)*logD-(k-1)*math.log(y+1))/((y*k*D/(D-1)+k)*logD)

				if less<gamma and more>gamma:

					r0 = y*k1

					possible+=[[D, r0]]

					possibleD+=[D]

					possibler0+=[r0]

				elif len(possible)>0 and possible[len(possible)-1][1]==y*k1:

					break

				D+=stepD*(Dmax-1)

			y+=stepy*y0

		for value in window.children.values():

			value.grid_remove()

		canvas = tkinter.Canvas(window, width = 480, height = 480, bg = 'white')

		canvas.grid(row = 0, column = 0)

		canvas.create_text(240, 470, text = 'Starting concentration')

		canvas.create_text(10, 240, angle = 90, text = 'Dilution coefficient (logarithmic scale)')

		canvas.create_line(40, 40, 40, 440)

		canvas.create_line(40, 440, 440, 440)

		f = open('bacteria_data.txt', 'w')

		canvas.create_text(30, 440, text = '1', angle = 90)

		startR1 = str(y0*k1//1)

		startR2 = str(y0*k1)[0:6]

		if len(startR1)>len(startR2):

			startR = startR1

		else:

			startR = startR2

		maxD1 = str(max(possibleD)//1)

		maxD2 = str(max(possibleD))[0:6]

		if len(maxD1)>len(maxD2):

			maxD = maxD1

		else:

			maxD = maxD2

		maxR1 = str(max(possibler0)//1)

		maxR2 = str(max(possibler0))[0:6]

		if len(maxR1)>len(maxR2):

			maxR = maxR1

		else:

			maxR = maxR2

		canvas.create_text(40, 450, text = startR)

		canvas.create_text(30, 40, text = maxD, angle = 90)

		canvas.create_text(440, 450, text = maxR)

		for dataset in possible:

			Daxis = 400*(math.log(dataset[0]))/(math.log(max(possibleD)))

			Raxis = 400*(dataset[1]-y0*k1)/(max(possibler0)-y0*k1)

			canvas.create_oval(39+Raxis, 439-Daxis, Raxis+41, 441-Daxis, outline = 'black', fill = 'black')

			f.write(str(dataset[0])+', '+str(dataset[1])+'\n')

		f.close()

		labelText = 'This plot is an approximation.\n'

		labelText+= 'Possible combinations written to file bacteria_data.txt.\n'

		labelText+= 'Every line of that file consists of two numbers\n'

		labelText+= 'separated by a comma and a space, where\n'

		labelText+= 'the first number is dilution coefficient,\n'

		labelText+= 'the second number is the corresponding\n'

		labelText+= 'starting concentration.\n'

		labelText+= 'The plot is clickable; if you click on a point in\n'

		labelText+= 'the plot area, corresponding points will be marked\n'

		labelText+= 'and labeled on the axes.'

		label = tkinter.Label(window, text = labelText)

		label.grid(row = 1, column = 0)

		button = tkinter.Button(window, text = 'Exit', command = destruction)

		button.grid(row = 2, column = 0)

		def canvasClick(event):

			if 40<=event.x and 440>=event.x and 40<=event.y and 440>=event.y:

				canvas.delete('text')

				canvas.create_line (event.x, 445, event.x, 435, tag = 'text')

				canvas.create_line (35, event.y, 45, event.y, tag = 'text')

				DaxisClick = 440-event.y 

				RaxisClick = event.x-40

				Raxis = RaxisClick*(max(possibler0)-y0*k1)/400+y0*k1

				Daxis = DaxisClick*math.log(max(possibleD))/400

				D = math.e**Daxis

				D1 = str(D//1)

				D2 = str(D)[0:6]

				if len(D1)>len(D2):

					textD = D1

				else:

					textD = D2

				canvas.create_text(30, event.y, text = textD, angle = 90, tag = 'text')

				R1 = str(Raxis//1)

				R2 = str(Raxis)[0:6]

				if len(R1)>len(R2):

					textR = R1

				else:

					textR = R2

				canvas.create_text(event.x, 450, text = textR, tag = 'text')

			return()

		canvas.bind('<Button-1>', canvasClick)

	return()

​

#Here you can change the in-program parameters. stepD is how much greater next dilution coefficient checked is than the previous, measured in fractions of Dmax-1

#(with dilution coefficient greater than Dmax coexistence is surely impossible)

stepD = 0.01

#stepy is how much greater next dimensionless resource concentration value is than the previous, measured in fractions of y0 

#(with value of y0 of less coexistence is surely impossible).

stepy = 0.1

#In line 37 you can appoint a maximum value for dimensionless resource concentration value y.

​

window = tkinter.Tk()

label = tkinter.Label(window, text = 'Enter the characteristics of your bacteria species.')

label.grid(row = 0, column = 0, columnspan = 2)

label1 = tkinter.Label(window, text = 'Species 1')

label1.grid(row = 1, column = 0, columnspan = 2)

labelg1 = tkinter.Label(window, text = 'Growth rate (g\u2081): ')

labelg1.grid(row = 2, column = 0)

entryg1 = tkinter.Entry(window)

entryg1.grid(row = 2, column = 1)

labelk1 = tkinter.Label(window, text = 'Monod constant (K\u2081): ')

labelk1.grid(row = 3, column = 0)

entryk1 = tkinter.Entry(window)

entryk1.grid(row = 3, column = 1)

label2 = tkinter.Label(window, text = 'Species 2')

label2.grid(row = 4, column = 0, columnspan = 2)

labelg2 = tkinter.Label(window, text = 'Growth rate (g\u2082): ')

labelg2.grid(row = 5, column = 0)

entryg2 = tkinter.Entry(window)

entryg2.grid(row = 5, column = 1)

labelk2 = tkinter.Label(window, text = 'Monod constant (K\u2082): ')

labelk2.grid(row = 6, column = 0)

entryk2 = tkinter.Entry(window)

entryk2.grid(row = 6, column = 1)

enterButton = tkinter.Button(window, text = 'OK', command = checkData)

enterButton.grid(row = 7, column = 0)

exitButton = tkinter.Button(window, text = 'Exit', command = destruction)

exitButton.grid(row = 7, column = 1)

window.mainloop()

