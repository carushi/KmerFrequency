#-*- coding:utf-8 -*-

import numpy as np

class BaseAnn(object):
	"""Calculate GC contents and Kmer frequency"""	
	"""A sequence is assumed to be encoded by 4 types of DNA base and 8 types of genome annotation"""
	"""and expressed by characters of seq"""
	"""(pos%4 -> ACGT, pos/4 -> each annotation)"""
	"""Ignore Ns"""	

	seq = "ACGTacgtopqrOPQRWXYZDEFHIJKLijkl"

	def __init__(self):
		self.rev = False
		pass

	def __del__(self):
		pass

	def setRev():
		self.rev = True
	def unsetRev():
		self.rev = False

	def getAnn(self, c):
		return self.seq.find(c)/4

	def getBase(self, c):
		return self.seq.find(c)%4

	def IsGC(self, c):
		n = self.getBase(c)
		return n == 1 or n == 2

	def CompInd(self, c):
		if c == 0:	return 3
		elif c == 1:	return 2
		elif c == 2:	return 1
		else:	return 0

	def kmerHash(self, fragment):
		list = []
		for c in fragment:
			list.append(self.getBase(c))
		num = 0
		if self.rev:
			for i in range(len(list)):
				if list[i] < 0: return -1
				num = num*4+self.CompInd(list[i])
		else:
			for i in range(len(list)):
				if list[i] < 0: return -1
				num = num*4+list[i]
		return num

	def getBaseKmer(self, fragment, kmer, j, num):
		if fragment[j] == "N":
			return -1
		elif num < 0:
			num = self.kmerHash(fragment[j:(j+kmer)])
		else:
			num = (num-self.getBase(fragment[j-1])*4**(kmer-1))*4+self.getBase(fragment[j+kmer-1])
		return num

	def setKmers(self, fragment, kmer, length = -1):
		if length < 0:	length = 4**kmer+1
		assert(length <= 2 or length == 4**kmer or length == 4**kmer+1)
		kmers = np.array([0.0]*(length))
		if length <= 2:
			GC = 0
			kmers[0] = len(filter(lambda x: self.IsGC(x), list(fragment)))
		else:
			num = -1
			for j in range(0, len(fragment)-kmer+1):
				num = self.getBaseKmer(fragment, kmer, j, num)
				if num >= 0:
					kmers[num] = kmers[num]+1
		if length == 2 or length == 4**kmer+1:	kmers[-1] = 1
		return kmers

	def getStrPos(self, fragment):
		list = [x for x in [self.getAnn(c) for c in fragment] if x >= 0]
		return np.argmax(np.bincount(list))

	def getStrPosName(self, fragment):
		num = self.getStrPos(fragment)
		return ["Intergenic", "Repeat", "Intron", "Coding Exon", "Non-coding exon", "Non-coding intron", "5'UTR", "3'UTR"][num]


		
if __name__ == '__main__':
	base = BaseAnn()
	print "4mer frequency"
	print base.setKmers("ACGTAAAACGTacgtopoooooppqqqrrrOPQRWXYZDEFHIJKLijkl", 4, 257)
	print "Major annotation : ",
	print base.getStrPosName("ACGTAAAACGTacgtopoooooppqqqrrrOPQRWXYZDEFHIJKLijkl")
