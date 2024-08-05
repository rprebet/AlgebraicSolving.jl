"""
25/02/2023
"""

import os
import subprocess

def PolyToBRstring(f):
	fc, fe = f.coefficients(), f.exponents()
	s = "{} 2\n".format(len(fc))
	M = [ "{} {} {}".format(fc[i], fe[i][1], fe[i][0]) for i in range(len(fc)) ]
	s += "\n".join(M)
	return s
	
def BRstringToPoly(s):
  Ls = s.split("\n") 
  Ls = [ ls.split(" ") for ls in Ls]
  Ls = [ [int(l) for l in ls] for ls in Ls ]
  
  A.<x,y> = ZZ[]
  F, i = [], 0
  while i<len(Ls):
  	  nterm = Ls[i][0]
  	  j, f = 1, A(0)
  	  while j<=nterm:
  	  	c = Ls[i+j]
  	  	f += c[0]*y^(c[1])*x^(c[2])
  	  	j+=1
  	  F.append(f)
  	  i += nterm + 1
  
  if len(F)==1:
  	return f
  else:
  	return F

def ToBresultant(f, g, finput="/tmp/in.bs"):
	"""Convert two bivariate polynomials into a bresultant input file.
		w.r.t. the first variable w.r.t. the ordering
	
	Inputs :
	f, g (polynomials): couple of bivariate polynomials we want to solve s.t x<y
	finput (string): path to the bresultant input file.

	"""	
	assert f.parent() == g.parent(), "The two polynomials must belong to the same polynomial ring."
	assert f.parent().characteristic() == 0, "Positive characteristic is not handled"

	s = PolyToBRstring(f) + "\n" + PolyToBRstring(g)
		
	fd = open(finput, 'w')
	fd.write(s)
	fd.close()
	
def FromBresultant(output="p", foutput="/tmp/out.bs"):
	"""Convert an output of bresultant to a sage object
	
	Inputs :
	output ["p", "r", "q"]: type of output
		* "p": parametrizations of the resultant sq-free factors
		* "r": resultant (todo)
		* "q": square-free part of the resultant (todo, increasing degree as output)

	"""	
	A.<x,y> = ZZ[]
	f = open(foutput,'r')
	s = f.read()
	s = s.replace("\n","").replace(":","").replace(";","")
	L = sage_eval(s, locals={'x': x, 'y': y})
	return L
	
def Bresultant(f, g, fname1="/tmp/in.bs", fname2="/tmp/out.bs", bspath="./bresultant", v=0):
	ToBresultant(f, g, fname1)
	out = subprocess.run([bspath, "-f", fname1, "-o", fname2], capture_output=True,text=True)
	if v==1:
		print((out.stderr).split('\n')[-2].split(': ')[-1]+'s')
	if v==2:
		print(out.stderr)
	return FromBresultant(fname2)
	
#A.<x,y> = ZZ[]                                                            
#f=A.random_element()                                                             
#g=diff(f,y)                                                      
#Bresultant(f,g,v=2)
#FromBresultant("p")
