# minitest.py
#
# Smart testing system.
#
# 2017 Alexander Oleynichenko
# alexvoleynichenko@gmail.com
#

#!/usr/bin/env python

import sys
import os
import re
import subprocess


# main options
VERBOSE = 0
FSCC_PATH = "minichem.x"


def clean():
	files = ["AOINTS*"]
	for f in files:
		os.system("rm -f " + f)

class Filter:
	
	def __init__(self, pattern, answer, eps):
		self.pat = pattern
		self.ans = [answer] if isinstance(answer, float) else answer
		self.eps = [eps] if isinstance(eps, float) else eps
		if len(self.eps) == 1 and len(self.ans) > 1:
			self.eps = self.eps * len(self.ans)
		self.re_float = re.compile(r'([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?')
		if len(self.eps) != len(self.ans):
			raise 'lendth of answers and epsilons should be equal'
	
	def __str__(self):
		if len(self.ans) == 1:
			return "%-22s%-20.8g%10.2e" % ('"'+self.pat+'"', self.ans[0], self.eps[0])
		else:
			return "%-22s%-20.8g%10.2e" % ('"'+self.pat+'"', self.ans[0], self.eps[0])
	
	def match(self, string):
		try:
			string = string.split(self.pat,1)[1]
			nums = re.findall(self.re_float, string)
			nums = [float(num[0]) for num in nums]
		except:
			return False
		if len(self.ans) > len(nums):
			return False
		
		# check for coincidence
		for i,e in enumerate(self.eps):
			if abs(nums[i]-self.ans[i]) >= e:
				return False
		return True


class Test:
	
	def __init__(self, name, inp, filters):
		self.inp = inp
		self.filters = filters
		self.name = name
	
	"""
	Returns True if OK, False if the test was failed.
	We are now in working directory with file 'name'
	"""
	def run(self):
		import time
		
		output_name = self.inp + ".test.out"
		matches = [False for f in self.filters]
		cmd = FSCC_PATH + " " + self.inp + " | tee " + output_name
		
		t1 = time.time()
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		out, err = p.communicate()
		t2 = time.time()
		result = out.split('\n')
		for lin in result:
			for i,f in enumerate(self.filters):
				if matches[i]:  # skip value which is already done
					continue
				matches[i] = f.match(lin)
		# remove temporary files
		clean()
		
		status = False not in matches
		print '    %-40s%10.3f sec %15s' % (self.name[:40], (t2-t1), 'PASSED' if status else 'FAILED')
		# print detailed info about failure
		if not status:
			print "\n      Failure diagnostics"
			print "      Input file : ", self.inp
			print "      Output file: ", output_name
			print "      !  Text pattern         Expected value        Epsilon          Status"
			for i,f in enumerate(self.filters):
				print "      ! %50s         %s" % (str(f), "Passed" if matches[i] else "Failed")
			print
		return status
