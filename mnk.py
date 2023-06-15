import numpy as np

def solve_mnk(x, y):
	"""
	returns MNK coofs in format (k, s_k, b, s_b)
	"""
	axx = (x*x).mean()
	ayy = (y*y).mean()
	axy = (x*y).mean()
	ax = x.mean()
	ay = y.mean()
	k = (axy-ax*ay)/(axx-ax*ax)
	b = ay-k*ax
	s_k = ((((ayy-ay*ay)/(axx-ax*ax))-k*k)/(len(x)-2))**0.5
	s_b = s_k*(axx**0.5)
	return (k, s_k, b, s_b)