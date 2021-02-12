import numpy as np
import math

class test_fx:
	def __init__(self, func, x_range):
		self.func = func
		self.x_range = x_range
		
		self.y_range = []
		self.slope_range = []
		
		self.num_roots = None
		self.root_intervals = None
		self.monotonic = None

		self.continuous = None
		self.num_singularities = None

		self.mean = None
		self.std = None

		self.dx = 1.e-6
		self.max_points = round(1.e7)
		self.interval_size, self.num_intervals, self.intervals = None, None, []

		self.update_sub_dividing()
		
		self.temp_interval_data = []

	def full_tests(self):
		return

	def update_sub_dividing(self, dx=self.dx):
		self.dx = dx

		self.interval_size = self.dx*self.max_points
		self.num_intervals = (self.x_range[1]-self.x_range[0])//self.interval_size
		self.num_intervals += 1 if self.num_intervals*self.interval_size < (self.x_range[1]-self.x_range[0]) else 0
		self.intervals = [[self.x_range[0]+i*self.interval_size, self.x_range[0]+(i+1)*self.interval_size] for i \
							in range(slef.num_intervals-1)]
		self.intervals += [[self.x_range[0]+(self.num_intervals-1)*self.interval_size, self.x_range[1]]]
		
		return

	def probe_around(self, x_0):
		count = 0
		max_attempts = 10
		x, y = None, None

		while count < max_attempts:
			try:
				x = x_0+(np.random.rand()-0.5)*self.dx 
				y = self.func(x)
			except:
				count += 1
			else:
				return (x, y)

		print("invalid function input")
		exit(1)

	def values(self, interval):
		N = round((interval[1]-interval[0])/self.dx) + 1
		value_pairs = []
		x, y = None, None

		for i in range(N):
			try:
				x = i*self.dx
				y = self.func(x)
			except:
				x, y = self.probe_around(i*self.dx)
			value_pairs.append([x, y])
		
		return value_pairs


	def monotonic(self, test_temp_interval=False):
		result = 1
		
		if test_temp_interval:
			for i in range(1, len(self.temp_interval_data)):








