#Python3 script to create random timetable instances
import random
import math
import pprint
import csv
import os, sys
from ortools.sat.python import cp_model
import itertools as it
import csv
from pandas import DataFrame, read_csv
import pandas as pd 
import copy
import pickle
import time
import networkx as nx
import matplotlib.pyplot as plt



def s2(x_m, x2_m, flist,threads,timeblock,h_mprev = None):
	print('s2')
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	
#	Sorig = [s for s in S if sd[s].course_name[4] in ('3')]	
	Sorig = list(S)

		
	SallLab = [s for s in Sorig if sd[s].isExtended if sd[s].periods > 1.5]	
	SallLab.append(0)
	

	S933 = [s for s in Sorig if sd[s].course_name[4] in ('4') if s not in SallLab]	
	SallLab.extend(S933)
	
	
	Sreg = [s for s in Sorig if s not in SallLab]
			
	if Lunch != -1:
		T = [t for t in T if t != Lunch ]
		
	Y = set([(s,d,t,r) for s in S for d in D for t in T for r in R if (sd[s].roomtype == rd[r].roomtype) ])
	W = set([(g,s) for g in G for s in S if (sd[s].course_num in gd[g].course_num_list) ]  )	
	Q = [(p,s) for p in P for s in S if sd[s].course_num in  pd[p].secNums ]
	
	S3 = [s for s in S if sd[s].isExtended and sd[s].periods == 3]
	
	if Lunch != -1:
		T3 = [t for t in T if ((t+2 < Lunch) or ((t > Lunch) and t+2 <= max(T) )) ]
	else:
		T3 = [t for t in T if t+2 <= max(T)  ]
		
	Y3 = set([(s,d,t,r) for s in S3 for d in D for t in T3 for r in R if (sd[s].roomtype == rd[r].roomtype) ])
	Y3h = set([(s,d,t) for s in S3 for d in D for t in T3 ])
	
	S4 = [s for s in S if sd[s].isExtended and sd[s].periods == 4]
	
	if Lunch != -1:
		T4 = [t for t in T if t+3 < Lunch ]
	else:
		T4 = [t for t in T if t+3 <= max(T) ]
		
	Y4 = set([(s,d,t,r) for s in S4 for d in D for t in T4 for r in R if (sd[s].roomtype == rd[r].roomtype) ])
	Y4h = set([(s,d,t) for s in S4 for d in D for t in T4 ])
	
	S2 = [s for s in S if sd[s].isExtended and sd[s].periods == 2]
	if Lunch != -1:
		T2 = [t for t in T if ((t+1 < Lunch) or ((t > Lunch) and t+1 <= max(T) )) ]
	else:
		T2 = [t for t in T if t+1 <= max(T) ]
		
	Y2 = set([(s,d,t,r) for s in S2 for d in D for t in T2 for r in R if (sd[s].roomtype == rd[r].roomtype) ])
	Y2h = set([(s,d,t) for s in S2 for d in D for t in T2 ])
	
	PS = [(p,s) for p in P for s in S if s in pd[p].secNums]
	PDT = [(p,d,t) for p in P for d in D for t in T]
	GDT = [(g,d,t) for g in G for d in D for t in T]
	
	PDTS = set([(p,d,t,s) for p in P for d in D for t in T for s in S if (p,s) in Q])
	GDTS = set([(g,d,t,s) for g in G for d in D for t in T for s in S if ( ( (g,s) in W)  )  ])
	GC = [(g,c) for g in G for c in C if c in gd[g].course_num_list ]	
	DTR = set([(d,t,r) for d in D for t in T for r in R])
	print('model started')
	m = cp_model.CpModel()
	SDT = [(s,d,t) for s in S for d in D for t in T]
	
	numRoomInType = {}
	for r in R:
		rt = rd[r].roomtype
		if rt in numRoomInType:
			numRoomInType[rt]+=1
		else:
			numRoomInType[rt]=1
	
	RT = list(numRoomInType.keys())			
	
	
	#################Testing
	SS = [(s,t) for s in S for t in S if s < t if x2_m[s,t] > 0.5]
	print('Num sections ',len(S))
	print('Num edges ', len(SS))

	
	
	solver = cp_model.CpSolver()
#	solver.parameters.log_search_progress = True
	solver.parameters.num_search_workers = int(threads)
	solver.parameters.max_time_in_seconds = timeblock
#	status = solver.Solve(m)
	solution_printer = cp_model.ObjectiveSolutionPrinter()
	m = cp_model.CpModel()


	h = {(s,d,t): m.NewBoolVar('h%i%i%i' % (s,d,t) ) for (s,d,t) in SDT}

	if h_mprev != None:
		for (s,d,t) in SDT:
			m.AddHint(h[s,d,t],h_mprev[s,d,t])
	
	#common period
	# (comd,comt) = (0,5)
# 
# 	m.Add(h[0,0,5] + h[0,0,3]  == 1)

	for s in S:
			m.Add(sum(h[s,d,t] for d in D for t in T if (s,d,t) in SDT ) == sd[s].periods)



	fe4h = {(s,d,t): m.NewBoolVar('fe4h%i%i%i' % (s,d,t) ) for (s,d,t) in Y4h}
	for s in S4:
		m.Add(sum( fe4h[s,d,t] for d in D for t in T if (s,d,t) in Y4h)  == 1 )


	for (s,d,t) in Y4h:
		m.Add(h[s,d,t] + h[s,d,t+1] + h[s,d,t+2] + h[s,d,t+3] >= 4*fe4h[s,d,t] )
	
	for (s,d,t) in Y4h:	
		m.AddImplication(fe4h[s,d,t],h[s,d,t] )
		m.AddImplication(fe4h[s,d,t],h[s,d,t+1] )
		m.AddImplication(fe4h[s,d,t],h[s,d,t+2] )
		m.AddImplication(fe4h[s,d,t],h[s,d,t+3] )


	fe3h = {(s,d,t): m.NewBoolVar('fe3h%i%i%i' % (s,d,t)) for (s,d,t) in Y3h } 

	for s in S3:
		m.Add(sum( fe3h[s,d,t] for d in D for t in T if (s,d,t) in Y3h)  == 1 )

	for (s,d,t) in Y3h:	
		m.AddImplication(fe3h[s,d,t],h[s,d,t] )
		m.AddImplication(fe3h[s,d,t],h[s,d,t+1] )
		m.AddImplication(fe3h[s,d,t],h[s,d,t+2] )

	fe2h = {(s,d,t): m.NewBoolVar('fe2h%i%i%i' % (s,d,t)) for (s,d,t) in Y2h }

	for s in S2:
		m.Add(sum( fe2h[s,d,t] for d in D for t in T if (s,d,t) in Y2h)  == 1 )

	for (s,d,t) in Y2h:	
		m.AddImplication(fe2h[s,d,t],h[s,d,t] )
		m.AddImplication(fe2h[s,d,t],h[s,d,t+1] )
	
	t2times = {(s,d):m.NewBoolVar('t2times%i%i' % (s,d) ) for s in S for d in D }
	
	for s in S: 
		for d in D:
			if sd[s].periods == 1 or not sd[s].isExtended:
				m.Add( sum( h[s,d,t] for t in T  if (s,d,t) in SDT ) <= 1 + t2times[s,d])	


	troomh = {(d,t,rt):m.NewBoolVar('troomh%i%i%s' % (d,t,rt) ) for d in D for t in T for rt in RT}
	for d in D:
		for t in T:
			for rt in RT:
				m.Add(sum( h[s,d,t] for s in S if sd[s].roomtype == rt) <= numRoomInType[rt] + 2*troomh[d,t,rt]) 


						
	
	SS = [(s,t) for s in S for t in S if s < t if x2_m[s,t] > 0.5]
	tconf = {(s1,s2,d,t):m.NewBoolVar('tconf%i%i%i%i' % (s1,s2,d,t) ) for (s1,s2) in SS for d in D for t in T }
	
	pof = {(p,d): m.NewBoolVar('pof%i%i' % (p,d)) for p in P for d in D}

	for p in P:
		m.Add( sum(pof[p,d] for d in D) <= len(D)-1 )
#		m.Add( pof[p,random.choice(D)] == 0 )
		m.Add( pof[p,p%5] == 0 )
	


	tpof = {p: m.NewBoolVar('tpof%i' % p) for p in P}
	for (p,s) in PS:
		for d in D:
			for t in T:
				m.Add(h[s,d,t] <= tpof[p]).OnlyEnforceIf(pof[p,d].Not())
	comw = {s: 1 for s in S}
	comw[0] = 10
	
	m.Minimize( 10*sum(t2times[s,d] for s in S for d in D)  + sum(tpof[p] for p in P ) + sum(1000*comw[s1]*tconf[s1,s2,d,t]  for (s1,s2) in SS for d in D for t in T) + 100*sum(troomh[d,t,rt] for d in D for t in T for rt in RT) )
	
	startI=0
	if h_mprev != None:
		startI = 0

	tabu1 = None
	for i in range(startI,2):
		if i == 0:
			limitTo = SallLab
			solver.parameters.max_time_in_seconds = timeblock/4	
			
		if i == 1:
			limitTo = Sorig	
			solver.parameters.max_time_in_seconds = (3/4)*timeblock
		for (s1,s2) in SS:
			for d in D:
				for t in T:
					if s1 in limitTo and s2 in limitTo:
						m.Add( h[s1,d,t] + h[s2,d,t] <= 1 + tconf[s1,s2,d,t] )
						
		print('Calling Solver')
		status = solver.SolveWithSolutionCallback(m, solution_printer)
		if status == cp_model.INFEASIBLE:
			print('Infeasible value of i is ',i)
			return None
		if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
			print('finished sections and common hour')
			Temp104 =  [(s1,s2,d,t) for (s1,s2) in SS for d in D for t in T if solver.Value(tconf[s1,s2,d,t]) > 0.5]
			if len(Temp104) > 0:
				tabu1 = []
				for (s1,s2,d,t) in Temp104:
					print('Conflict: ',s1, s2, sd[s1].course_name, sd[s2].course_name, ' time date ',d,t)
					#SS = [(s,t) for s in S for t in S if s < t if x2_m[s,t] > 0.5]
					Temp105 = [g for g in G if (g,s1) in W if (g,s2) in W if x_m[g,s1] > 0.5 if x_m[g,s2] > 0.5]
					for g in Temp105:
						tabu1.append((g,s1))
						tabu1.append((g,s2))
					
					
			U1Ah = [(s,d,t) for (s,d,t) in SDT if s in limitTo if solver.Value(h[s,d,t]) > 0.5]
			U1Bh = [(s,d,t) for (s,d,t) in SDT if s in limitTo if solver.Value(h[s,d,t]) < 0.5]
			
			for q in U1Ah:
				m.AddHint(h[q],1)
			for q in U1Bh:
				m.AddHint(h[q],0)
			
			if i==1:	
				for q in U1Ah:
					m.Add(h[q]==1)
				for q in U1Bh:
					m.Add(h[q]==0)		
			
	
		
			


	
#	status = solver.SolveWithSolutionCallback(m, solution_printer)
#	solver.parameters.max_time_in_seconds = 100000

	if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
		z = {(s,d,t,r): m.NewBoolVar('z%i%i%i%i' % (s,d,t,r) ) for (s,d,t,r) in Y}
		for (s,d,t) in SDT:					
			m.Add( h[s,d,t] == sum(z[s,d,t,r] for r in R if (s,d,t,r) in Y)  )
		for s in S:
			m.Add(sum(z[s,d,t,r] for d in D for t in T for r in R if (s,d,t,r) in Y ) == sd[s].periods)
		
		fe4 = {(s,d,t,r): m.NewBoolVar('fe4%i%i%i%i' % (s,d,t,r) ) for (s,d,t,r) in Y4}
		for s in S4:
			m.Add(sum( fe4[s,d,t,r] for (d,t,r) in DTR if (s,d,t,r) in Y4)  == 1 )
		
		for (s,d,t,r) in Y4:
			m.Add(z[s,d,t,r] + z[s,d,t+1,r] + z[s,d,t+2,r] + z[s,d,t+3,r] >= 4*fe4[s,d,t,r] )
		
		for (s,d,t,r) in Y4:	
			m.AddImplication(fe4[s,d,t,r],z[s,d,t,r] )
			m.AddImplication(fe4[s,d,t,r],z[s,d,t+1,r] )
			m.AddImplication(fe4[s,d,t,r],z[s,d,t+2,r] )
			m.AddImplication(fe4[s,d,t,r],z[s,d,t+3,r] )
		
		
		fe3 = {(s,d,t,r): m.NewBoolVar('fe3%i%i%i%i' % (s,d,t,r)) for (s,d,t,r) in Y3 } 
		
		for s in S3:
			m.Add(sum( fe3[s,d,t,r] for (d,t,r) in DTR if (s,d,t,r) in Y3)  == 1 )

		for (s,d,t,r) in Y3:	
			m.AddImplication(fe3[s,d,t,r],z[s,d,t,r] )
			m.AddImplication(fe3[s,d,t,r],z[s,d,t+1,r] )
			m.AddImplication(fe3[s,d,t,r],z[s,d,t+2,r] )
		
		fe2 = {(s,d,t,r): m.NewBoolVar('fe2%i%i%i%i' % (s,d,t,r)) for (s,d,t,r) in Y2 }

	
		for s in S2:
			m.Add(sum( fe2[s,d,t,r] for (d,t,r) in DTR if (s,d,t,r) in Y2)  == 1 )

		for (s,d,t,r) in Y2:	
			m.AddImplication(fe2[s,d,t,r],z[s,d,t,r] )
			m.AddImplication(fe2[s,d,t,r],z[s,d,t+1,r] )
		
		
		troom = {(d,t,r):m.NewBoolVar('troom%i%i%i' % (d,t,r) ) for d in D for t in T for r in R}
		for d in D:
			for t in T:
				for r in R:
					m.Add(sum( z[s,d,t,r] for s in S if (s,d,t,r) in Y ) <= 1 + troom[d,t,r])
		
		
		m.Minimize( 10*sum(t2times[s,d] for s in S for d in D)  + sum(tpof[p] for p in P ) + sum(1000*tconf[s1,s2,d,t]  for (s1,s2) in SS for d in D for t in T) + 100*sum(troom[d,t,r] for d in D for t in T for r in R) )
		
		print('Calling Solver')
		status = solver.SolveWithSolutionCallback(m, solution_printer)
		if status == cp_model.INFEASIBLE:
			print('Infeasible Rooms')
			return
		if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
			print('reading z')
		
			print('reading z')
			z_m = {}
			for q in Y:
				z_m[q] = solver.Value(z[q])
	# 			if z_m[q] > 0.5:
	# 				print(q,z_m[q])
			print('computing w')

			h_m = {(s,d,t): solver.Value(h[s,d,t]) for s in S for d in D for t in T }
			w_m = {}
			for (p,d,t) in PDT:
				w_m[p,d,t] = sum([h_m[s,d,t] for s in S if (p,s) in PS ] )
	
	# 		
	# 		for (p,d,t) in PDT:
	# 			w_m[p,d,t] = sum([z_m[s,d,t,r] for s in S for r in R if (p,s) in PS if (s,d,t,r) in Y]) 
			print('computing u')
			u_m = {(g,d,t,s) : 0 for (g,d,t,s) in GDTS}


			for (g,s) in W:
				if x_m[g,s] > 0.5:
					for d in D:
						for t in T:
							u_m[g,d,t,s] = sum( [ z_m[s,d,t,r] for r in R if (s,d,t,r) in Y ] ) 
		
			# for (g,d,t,s) in GDTS:
	# 			u_m[g,d,t,s] = x_m[g,s]*sum( [ z_m[s,d,t,r] for r in R if (s,d,t,r) in Y ] ) 
	# 		
			print('creating pdf files')
		
			savedSol = (flist,x_m,z_m,w_m,u_m,GDTS)
			pickle.dump(savedSol, open("solX.p","wb") ) 
			return (h_m,tabu1,solver.ObjectiveValue())
 			
		else:
			print('Did not solve')		











def CreateMiniGroups(x_m,G,W,S):
	gr = 0
	mini1 = {}
	mini2 = {}
	mini3 = {}
	for g in G:
		c = set([s for s in S if (g,s) in W if x_m[g,s] > 0.5])	
		found = False
		for u in range(gr):
			if c == mini2[u]:
				mini1[u]+=1
				found = True
				mini3[u].append(g)
				break
		if not found:
			mini2[gr] = c
			mini1[gr] = 1
			mini3[gr] = [g]
			gr+=1
	
	print('Number of groups is ',gr)	
		
	return (gr,mini1,mini2,mini3)			
					
				




def greedy(flist):
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	GG = [(g1,g2) for g1 in G for g2 in G]
	for g in G:
		gd[g].course_num_set = set(gd[g].course_num_list)
		
	dmat = {(g,h): len((gd[g].course_num_set).symmetric_difference((gd[h].course_num_set))) for (g,h) in GG }
	
	secOfCourse = {}
	for c in C:
		secOfCourse[c] = [s for s in S if sd[s].course_num == c]
	
	
	enrolled = {}
	rostor = {}
	restOfFam = {}
	
	for s in S:
		enrolled[s] = 0
		rostor[s] = []
		c = sd[s].course_num
		if catalog[c].partOfFam:
			restOfFam[s] = [r for r in S if sd[s].labtie == sd[r].labtie if r != s ]
	
	
	#set up conflict graph
	E = set([])	
	SS = [(s1,s2) for s1 in S for s2 in S if s1 < s2]
	PS = [(p,s) for p in P for s in S if s in pd[p].secNums]
	PSS = [(p,s1,s2) for p in P for (s1,s2) in SS if ( (p,s1) in PS and (p,s2) in PS )  ]
	for (p,s1,s2) in PSS:
		E.add((s1,s2))
		E.add((s2,s1))	
	
	rtCount = {}
	for r in R:
		if rd[r].roomtype in rtCount:
			rtCount[rd[r].roomtype]+=1
		else:
			rtCount[rd[r].roomtype] = 1
			
		
	
	SoneRoom = [s for s in S if rtCount[sd[s].roomtype] == 1 ]
	
	RSS = [(r,s1,s2) for r in R for s1 in SoneRoom for s2 in SoneRoom if s1 < s2 if sd[s1].roomtype == sd[s2].roomtype]
	
	for (r,s1,s2) in RSS:
		E.add((s1,s2))
		E.add((s2,s1))	
		
	Ecop = copy.deepcopy(E)
	
		
	random.shuffle(G)
	
	gdcopy = copy.deepcopy(gd)
	
	
	
	enrolledcopy = copy.deepcopy(enrolled)
	rostorcopy = copy.deepcopy(rostor)
	
	
	for i in range(1):
		x = {}
		gsch = {}
		random.shuffle(G)
		E = copy.deepcopy(Ecop)
		gd = copy.deepcopy(gdcopy)
		enrolled = copy.deepcopy(enrolledcopy)
		rostor = copy.deepcopy(rostorcopy)
		#enroll the first student
		g = G[0]
		gsch[g] = []
		while len(gd[g].course_num_set) > 0.5:
			c = random.choice(tuple(gd[g].course_num_set))
			s = secOfCourse[c][0]
			x[g,s] = 1
			enrolled[s]+=1
			rostor[s].append(g)
			gsch[g].append(s)
			gd[g].course_num_set.remove(c)
			if catalog[c].partOfFam:
				for r in restOfFam[s]:
					u = sd[r].course_num
					x[g,r] = 1
					enrolled[r]+=1
					rostor[r].append(g)
					gsch[g].append(r)
					gd[g].course_num_set.remove(u)
	
		Eg = [(a,b) for a in gsch[g] for b in gsch[g] if a != b  ]
		E.update(Eg)		
		H = set(G)
		H.remove(g)
		while (len(H) > 0.5):
			nDist = [(dmat[(g,g2)],g2) for g2 in tuple(H)]
			nDist.sort(key=lambda x:x[0])
			h = nDist[0][1]
			gsch[h] = []
			T1 = [s for s in gsch[g] if sd[s].course_num in gd[h].course_num_set if enrolled[s] < sd[s].cap ]
			for s in T1:
				x[h,s] = 1
				enrolled[s]+=1
				rostor[s].append(h)
				gsch[h].append(s)
				c = sd[s].course_num
				gd[h].course_num_set.remove(c)
			while len(gd[h].course_num_set) > 0.5:	
				c = random.choice(tuple(gd[h].course_num_set))
				T1 = [s for s in secOfCourse[c] if enrolled[s] < sd[s].cap ]
				pp = 100000
				best = -1
				for s in T1:
					T2 = list(gsch[h])
					T2.append(s)
					Ehs = set([(a,b) for a in T2 for b in T2 if a != b  ])
					Rh = Ehs - E
					if len(Rh) < pp:
						pp = len(Rh)
						best = s
				s = best
				x[h,s] = 1
				enrolled[s]+=1
				rostor[s].append(h)
				gsch[h].append(s)
				gd[h].course_num_set.remove(c)
				if catalog[c].partOfFam:
					for r in restOfFam[s]:
						u = sd[r].course_num
						x[h,r] = 1
						enrolled[r]+=1
						rostor[r].append(h)
						gsch[h].append(r)
						gd[h].course_num_set.remove(u)
	
			
						
			Eh = [(a,b) for a in gsch[h] for b in gsch[h] if a != b  ]
			E.update(Eh)		
			H.remove(h)		
			g = h
		
		
		if i==0:
			best1 = len(E)/2
			xbest = x
			Ebest = E
		if i > 0:
			if len(E)/2 < best1:
				best1 = len(E)/2
				xbest = x	
				Ebest = E
	
	x2 = {}
	for (s1,s2) in SS:
		if (s1,s2) in Ebest:
			x2[s1,s2] = 1
		else:
			x2[s1,s2] = 0
	print(best1)
	return (xbest,x2)	
	
		



def secCP(threads, flist, time,x_m = None,x2_m = None, tabu1 = None):
	print('USMMA Sectionize')
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
#	G = G[:40]
	W = set([(g,s) for g in G for s in S if (sd[s].course_num in gd[g].course_num_list) ]  )
	GC = [(g,c) for g in G for c in C if c in gd[g].course_num_list ]	
#	Q = [(p,s) for p in P for s in S if sd[s].course_num in  pd[p].secNums ]
	SS = [(s1,s2) for s1 in S for s2 in S if s1 < s2]
	GSS = [(g,s1,s2) for g in G for (s1,s2) in SS if ( (g,s1) in W and (g,s2) in W )  ]
	PS = [(p,s) for p in P for s in S if s in pd[p].secNums]
	PSS = [(p,s1,s2) for p in P for (s1,s2) in SS if ( (p,s1) in PS and (p,s2) in PS )  ]
	labtieSS = [(s1,s2) for s1 in S for s2 in S if sd[s1].labtie != '' and sd[s1].labtie == sd[s2].labtie and sd[s1].iAmParent and s1 != s2 ]
	GlabtieSS = [(g,s1,s2) for g in G for (s1,s2) in labtieSS if (g,s1) in W]
	
	
	# for g in G:
# 		print(g,gd[g].course_num_list)
# 	
# 	
	C = set([c for (g,c) in GC])
	
	C = list(C)
	C.sort()
	
# 	for c in C:
# 		print(c,catalog[c].name, sum([sd[s].cap for s in S if sd[s].course_num == c]) - sum([1 for g in G if (g,c) in GC]))
# 	


	
	
	
	
	model = cp_model.CpModel()
	x = {(g,s) : model.NewBoolVar('x%i%i' % (g,s) ) for (g,s) in W}
	# x = {}
# 	for g in G:
# 		for s in S:
# 			if (g,s) in W:
# 				x[g,s] = model.NewBoolVar('x%i%i' % (g,s) )
	
	print('Building Model')
	x2 = {}
# 	for s in S:
# 		for t in S:
# 			if (s,t) in SS:
# 				x2[s,t] = model.NewBoolVar('x%i%i' % (s,t) )
			
	for (s,t) in SS:
		x2[s,t] = model.NewBoolVar('x%i%i' % (s,t) )
	
# 	x2 = [model.NewBoolVar('x%i%i' % (s,t) ) for (s,t) in SS]		
	
	weight1 = {}
	for (s1,s2) in SS:
		wei = 1
		if sd[s1].isExtended:
			wei+=3
		if sd[s2].isExtended:
			wei+=3	
		weight1[s1,s2] = wei
	
	if tabu1 != None:
		model.Minimize( sum(weight1[s,t]*x2[s,t] for (s,t) in SS ) + sum(5*x[g,s] for (g,s) in tabu1))
	else:
		model.Minimize( sum(x2[s,t] for (s,t) in SS ) )
	
	if x_m != None:
		for (g,s) in W:
			if (g,s) in x_m.keys():
				if x_m[g,s] == 1:
					model.AddHint(x[g,s],1)
			
		for (g,s) in W:
			if (g,s) not in x_m.keys():	
				model.AddHint(x[g,s],0)
		
		for (s1,s2) in SS:
			model.AddHint(x2[s1,s2], x2_m[s1,s2])
		
				
#			model.Add( x[g,s] == 1)
	
	rtCount = {}
	for r in R:
		if rd[r].roomtype in rtCount:
			rtCount[rd[r].roomtype]+=1
		else:
			rtCount[rd[r].roomtype] = 1
			
		
	
	SoneRoom = [s for s in S if rtCount[sd[s].roomtype] == 1 ]
	# for s in SoneRoom:
	# 	print(s,sd[s].course_name)
	
	RSS = [(r,s1,s2) for r in R for s1 in SoneRoom for s2 in SoneRoom if s1 < s2 if sd[s1].roomtype == sd[s2].roomtype]
	
	for (r,s1,s2) in RSS:
		model.Add( x2[s1,s2] == 1 )
	
	for (p,s1,s2) in PSS:
		model.Add( x2[s1,s2] == 1 )
# 		
	for (g,c) in GC:
		model.Add(sum(x[g,s] for s in S if sd[s].course_num == c ) == 1  )
# 	
	for s in S:
		model.Add(sum(x[g,s] for g in G if (g,s) in W ) <= sd[s].cap)	
# 	
	
	



	
	for (g,s1,s2) in GlabtieSS:	
		model.Add(x[g,s2] <= x[g,s1])
		
	for (g,s1,s2) in GSS:
		model.AddBoolOr([x[g,s1].Not(), x[g,s2].Not(), x2[s1,s2]])		
#		model.Add(x[g,s1] + x[g,s2] - 1 <= x2[s1,s2])
	
	
	
	
	print('Ready to Solve')
	solver = cp_model.CpSolver()
	print('Created Solver')
#	solver.parameters.log_search_progress = True

	solver.parameters.search_branching = (cp_model.sat_parameters_pb2.SatParameters.PORTFOLIO_WITH_QUICK_RESTART_SEARCH)
	
	solver.parameters.linearization_level = 0

	solver.parameters.max_time_in_seconds = time
	solver.parameters.num_search_workers = int(threads)
#	status = solver.Solve(model)
	solution_printer = cp_model.ObjectiveSolutionPrinter()
	status = solver.SolveWithSolutionCallback(model, solution_printer)
	print('Finished Solving')
	if status == cp_model.OPTIMAL:
		print('found optimal')
		print('Total weight = %i' % solver.ObjectiveValue())
    
	if status == cp_model.INFEASIBLE:
		print('Infeasible')	
		exit()
     
	if status == cp_model.FEASIBLE:
		print('found feasible')	
		print('Total weight = %i' % solver.ObjectiveValue())

	if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
		x_m = {}
		x2_m = {}
		for (s,t) in SS:
			x2_m[s,t] = solver.Value(x2[s,t])
		
		print('***Number of Edges ', len([(s,t) for (s,t) in SS if x2_m[s,t] > 0.5]))
		
		
		for (g,s) in W:
			x_m[g,s] = solver.Value(x[g,s])
			
		print('saving pickle')
		mytuple = (x_m,x2_m,flist)
		pickle.dump(mytuple, open( "saveX.p", "wb"))
		
		
		(gr,mini1,mini2,mini3) = CreateMiniGroups(x_m,G,W,S)
		return mytuple 	
			
			
	


def prepInstanceForSolvers(plist):
	(divDict,term2Requests,divSizeDict,catalog,courseNametoNum,rd,roomTypeList,roomTypeNameToNum,gd,sd,pd) = plist
	S = list(sd.keys())
	S.sort()
	R = list(rd.keys())
	R.sort()
	RT = [roomTypeNameToNum[rtn] for rtn in roomTypeList]
	RT.sort()
	C = list(catalog.keys())
	C.sort()
	P = list(pd.keys())
	P.sort()
	G = list(gd.keys())
	G.sort()
	D = [0,1,2,3,4]
	T = [0,1,2,3,4,5,6,7]
	Lunch = 4
	
	Labties = set([sd[s].labtie for s in S])
	L = list(Labties)
	L.sort()
	
	
	for s in S:
		c = sd[s].course_num
		catalog[c].secList.append(s)
	
	
	
	plist = (S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict)
	return plist
	
	



def makeProfs(plist):
	(divDict,term2Requests,divSizeDict,catalog,courseNametoNum,roomDict,roomTypeList,roomTypeNameToNum,gdict,sdict) = plist
	S = list(sdict.keys())
	S.sort()
	
	#figure out disciplines
	for s in S:
		v = sdict[s]
		name = v.course_name
		discp = name[:4]
		v.discp = discp
		v.load = v.periods
		if v.isExtended and v.periods > 1.5:
			v.load = v.periods*(0.66)
		
		
	pd = {}
	profNum = 0
	DIS = set([sdict[s].discp for s in S ])
	for d in DIS:
		Sd = set([s for s in S if sdict[s].discp == d])
		while len(Sd) > 0.5:
			s = Sd.pop()
			suitProf = [p for p in pd.keys() if (pd[p].discp == d) and (pd[p].load + sdict[s].load <= 12)  ]
			if len(suitProf) > 0.5:
				p = suitProf[0]
				pd[p].addSec(sdict[s])
			else:
				pd[profNum] = Prof(profNum,d)
				pd[profNum].addSec(sdict[s])
				profNum+=1
	
	return pd
	
				
				
				
	
		
		






def readCSVPandas(curricu):
	location = curricu + '/2024Cur/'
#	location = r'curriculumH/2024Cur/'
	courses = pd.read_csv(location + 'COURSES-Table 1.csv')
	currics = pd.read_csv(location + 'CURRICULUM-Table 1.csv')
	rooms =  pd.read_csv(location + 'ROOMS-Table 1.csv')
	divsizes = pd.read_csv(location + 'DIVSIZES-Table 1.csv')

	c = list(currics.columns.values)
	c.pop(0) # we don't need the first two columns
	c.pop(0)
	#print(c)
	CourseRequests = set([])
	for i in c:
		xcol = list(currics[i])
		for n in xcol:
			if pd.notna(n):
				CourseRequests.add(n.strip().upper())

			
	courseCat = set([])
	for n in courses["COURSE"]:
		if pd.notna(n):
			courseCat.add(n.strip().upper())

	missingCourses = CourseRequests - courseCat

	if len(missingCourses) > 0:
		print(missingCourses)
		exit('Requested courses are not in the course catalogue')
	roomNameToNum = {}
	roomTypeNameToNum = {}
	CC = list(courses.index)
	RR = list(rooms.index)
	roomTypeSet = set([])
	roomDict = {}
	for row in RR:
		name = rooms.loc[row,"ROOMNAME"].strip().upper()
		roomtype = rooms.loc[row,"GENTYPE"].strip().upper()
		roomTypeSet.add(roomtype) 	
		cap = rooms.loc[row,"ROOMCAP"]
		roomDict[row] = Room(name,row,roomtype,cap)
		roomNameToNum[name] = row
	
	
	roomTypeList = list(roomTypeSet)
	for i in range(len(roomTypeList)):
		roomTypeNameToNum[roomTypeList[i]] = i
		
	courseNametoNum = {}
	catalog = {}
	for row in CC:
		if pd.notna(courses.loc[row,"COURSE"]):
			name = courses.loc[row,"COURSE"].strip().upper()
			periods = courses.loc[row,"PERIODS"]
			roomtype = courses.loc[row,"ROOMTYPE"].strip().upper()
			if roomtype not in roomTypeSet:
				print(roomtype)
				exit('Roomtype not in database')
			cap = courses.loc[row,"CAP"]
			isExtended = False
			
			if pd.notna(courses.loc[row,"EXTENDED"]) and courses.loc[row,"EXTENDED"].strip() in ('y','Y'):
				isExtended = True
			if pd.notna(courses.loc[row,"PARENT"]):
				parent = courses.loc[row,"PARENT"].strip().upper()
				if parent not in courseCat:
					exit('Parent course not in catalogue')
			else:
				parent = ''
			catalog[row] = Course(name,row,isExtended,periods,roomtype,cap,parent)
			courseNametoNum[name] = row
	
			
	
	#get term 2 courses
	
	
	
	catDict = {}
	
	reqDict = {}
	
	DS = list(divsizes.index)
	divSizeDict = {}
	for row in DS:
		if currics.loc[row,"TERM"] == 2:
			divname = divsizes.loc[row,"DIVISION"].strip().upper()
			divSizeDict[divname] = divsizes.loc[row,"SIZE"]
			
			
			
		
	
	COLS = list(currics.columns.values)
	COLS.pop(0) # we don't need the first two columns
	COLS.pop(0)
	term2Requests = set([])
	
	divDict = {}
	
	CU = list(currics.index)

	
	
	divNum=0
	for r in CU:
		if currics.loc[r,"TERM"] == 2:
			mycurric = set([])
			numCurric = set([])
			divname = currics.loc[r,"DIVISION"].strip().upper()
			divsize = divSizeDict[divname]
			for c in COLS:
				cname = currics.loc[r,c]
				if pd.notna(cname):
					term2Requests.add(cname.strip().upper())
					mycurric.add(cname.strip().upper())
					numCurric.add(courseNametoNum[cname])
			divDict[divNum]  = Division(divname,divNum, currics.loc[r,"TERM"],divsize,list(mycurric), list(numCurric))	
			divNum+=1	
	
			
	plist = (divDict,term2Requests,divSizeDict,catalog,courseNametoNum,roomDict,roomTypeList,roomTypeNameToNum)				
	return plist

def createStudentsAndSections(plist):		
	(divDict,term2Requests,divSizeDict,catalog,courseNametoNum,roomDict,roomTypeList,roomTypeNameToNum)	= plist
	gdict = {}
	sdict = {}
	
	numStud = 0
	for (k,v) in divDict.items():
		for n in range(v.size):
			gdict[numStud] = Student(numStud,v.name+'.'+str(n),  v.name,v.number,v.course_name_list,v.course_num_list)
			gdict[numStud].numPeriods = sum([catalog[c].periods for c in v.course_num_list])
			# if gdict[numStud].numPeriods > 29:
			# 	print(v.name,gdict[numStud].numPeriods)
			numStud+=1
	
	

	
	#figure out num of sections needed
	freqDict = {}
	for (k,v) in gdict.items():
		for n in v.course_num_list:
			if n in freqDict:
				freqDict[n]+=1
			else:
				freqDict[n]=1
	
	C = list(freqDict.keys())
	C.sort()
	
	families0 = []
	Parents = set([catalog[c].parent for c in C if catalog[c].parent!='' ])
	for p in Parents:
		n = courseNametoNum[p]
		catalog[n].iAmParent = True
		b = [n] + [catalog[c].number for c in C if catalog[c].parent == p]
#		a = [p] + [(catalog[c].name,catalog[c].cap) for c in C if catalog[c].parent == p ]
#		n = courseNametoNum[p]
		families0.append(b)
	
	sats = []
	fams = []
	Cset = set(C)
#	print(families0)
	for f in families0:
		if max([catalog[c].cap for c in f]) != min([catalog[c].cap for c in f]):
			sats.append(f)
			for c4 in f:
				catalog[c4].partOfSat = True
			for c in f:
				Cset.remove(c)		
		else:
				fams.append(f)
				for c4 in f:
					catalog[c4].partOfFam = True
				
				
				for c in f:	
					Cset.remove(c)
	
	lecs = list(Cset)
	lecs.sort()		

			
	
	
#regular lectures without parents	
	iAmParent = False
	labtie = -1
	secNum = 0
	for c in lecs:
		v = catalog[c]
		numOfSec = math.ceil(freqDict[c]/v.cap)
		for n in range(numOfSec):
			newCap = math.ceil(freqDict[c]/numOfSec) + math.ceil(freqDict[c]/numOfSec)%2
#			newCap = math.ceil(freqDict[c]/numOfSec)
			sdict[secNum] = Section(secNum,v.isExtended,v.periods,v.roomtype,newCap,v.number,v.name,v.parent,v.name + '.' + str(n),labtie, iAmParent)
			secNum+=1
		
	for s in sats:
		labtie+=1
		n=0
		m = s[0]
		v = catalog[m]
		newCap = freqDict[m]
		iAmParent = True
		sdict[secNum] = Section(secNum,v.isExtended,v.periods,v.roomtype,newCap,v.number,v.name,v.parent,v.name + '.' + str(n),labtie, iAmParent)
		secNum+=1
		
		for i in range(1,len(s)):
			iAmParent = False
			c = s[i]
			v = catalog[c]
		numOfSec = math.ceil(freqDict[c]/v.cap)
		for n in range(numOfSec):
			newCap = math.ceil(freqDict[c]/numOfSec) + math.ceil(freqDict[c]/numOfSec)%2
			sdict[secNum] = Section(secNum,v.isExtended,v.periods, v.roomtype,newCap,v.number,v.name,v.parent,v.name + '.' + str(n),labtie, iAmParent)
			secNum+=1
			
			
	
	for f in fams:
		n=0
		c = f[0]
		iAmParent = True
		v = catalog[c]
		numOfSec = math.ceil(freqDict[c]/v.cap)
		newCap = math.ceil(freqDict[c]/numOfSec) + math.ceil(freqDict[c]/numOfSec)%2
		for n in range(numOfSec):
			labtie+=1
			for i in range(len(f)):
				c = f[i]
				v = catalog[c]
				if i == 0:
					iAmParent = True
				else:
					iAmParent = False	
				sdict[secNum] = Section(secNum,v.isExtended,v.periods, v.roomtype,newCap,v.number,v.name,v.parent,v.name + '.' + str(n),labtie, iAmParent)
				secNum+=1
	
	
	
	
	return (gdict,sdict) 	
	
	
	
	
	
	

def table_to_tuple(table):
	a = []
	for i in range(5):
		for j in range(8):
			a = a + [table[i][j]]
	return tuple(a)	


def createStats(flist,x,z,w,u,GDTS):
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	#enrollment per section 
	secEnroll = {}
	groupsInSec = {}
	for s in S:
		groupsInSec[s]= [g for g in G if (g,s) in x and x[g,s] > 0.5 ]
		secEnroll[s] = len( groupsInSec[s] )
	
	#room assigned to section, day, time
	roomAssign = {}
	for s in S:
		for d in D:
			for t in T:
				tmpR19 = [r for r in R if (s,d,t,r) in z and z[s,d,t,r] > 0.5 ]
				if len(tmpR19) > 1:
					exit('double scheduled rooms')
				if len(tmpR19) == 1:
					roomAssign[s,d,t] = tmpR19[0]
	
	return [secEnroll, roomAssign, groupsInSec]		


def makeRoomTables(flist,x,z,w,u,GDTS):
	[secEnroll, roomAssign, groupsInSec] = createStats(flist,x,z,w,u,GDTS)
	gtt = {}
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	
	list_of_colors = ['\\cellcolor{purple!30}', '\\cellcolor{blue!30}', '\\cellcolor{green!30}', '\\cellcolor{orange!30}', '\\cellcolor{gray!30}', '\\cellcolor{purple!60}', '\\cellcolor{blue!60}', '\\cellcolor{green!60}','\\cellcolor{orange!60}', '\\cellcolor{gray!60}', '\\cellcolor{purple!90}', '\\cellcolor{blue!90}', '\\cellcolor{green!90}','\\cellcolor{orange!90}', '\\cellcolor{gray!90}', '\\cellcolor{pink!30}','\\cellcolor{pink!60}', '\\cellcolor{pink!90}'  ]
	
	
	
# 	Q = [(p,s) for p in P for s in S if sd[s].course_name in  pd[p].course_numbers_list ]
# 	PS = [(p,s) for p in P for s in S if s in pd[p].secs]
	
	
	color = {}
	for s in S:
		color[s] = list_of_colors[s%len(list_of_colors)]
		
	
	
	for r in R:		
		table = [['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', '']]
		for d in D:
			for t in T:
				for s in S:
					if sd[s].course_name[2] == '&':
						sd[s].course_name = sd[s].course_name[:2] + sd[s].course_name[3:]	
					if ((s,d,t,r) in z) and (z[s,d,t,r] > 0.5):
						table[d][t] = table[d][t] + color[s] + 'Section ' + str(s) + " Extended "  +str(sd[s].isExtended) + " cap " + str(sd[s].cap) + " roomtype " + str(sd[s].roomtype) + " course " +  str(sd[s].course_name) + " enrol=" + str(str(secEnroll[s])) + " room=" + str(str(roomAssign[s,d,t]))
		gtt[r] = table
	return gtt	
		



				


def makeSecTables(plist,x,z,w,u,GDTS):
	[secEnroll, roomAssign, groupsInSec] = createStats(plist,x,z,w,u,GDTS)
	gtt = {}
	[S, sd,P,pd,G,gd,C,cd,roomtypeList, roomMaster, R, rd, roomToType, D,T, Lunch] = plist
	
	list_of_colors = ['\\cellcolor{purple!30}', '\\cellcolor{blue!30}', '\\cellcolor{green!30}', '\\cellcolor{orange!30}', '\\cellcolor{gray!30}', '\\cellcolor{purple!60}', '\\cellcolor{blue!60}', '\\cellcolor{green!60}','\\cellcolor{orange!60}', '\\cellcolor{gray!60}', '\\cellcolor{purple!90}', '\\cellcolor{blue!90}', '\\cellcolor{green!90}','\\cellcolor{orange!90}', '\\cellcolor{gray!90}', '\\cellcolor{pink!30}','\\cellcolor{pink!60}', '\\cellcolor{pink!90}'  ]
	
	
	
	Q = [(p,s) for p in P for s in S if sd[s].course_name in  pd[p].course_numbers_list ]
	PS = [(p,s) for p in P for s in S if s in pd[p].secs]
	
	
	for s in S:
		color = {}
		color[s] = list_of_colors[0]
				
		table = [['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', '']]
		for d in D:
			for t in T:
				for r in R:
					if ((s,d,t,r) in z) and (z[s,d,t,r] > 0.5):
						table[d][t] = table[d][t] + color[s] + 'Section ' + str(s) + " Extended "  +str(sd[s].isExtended) + " cap " + str(sd[s].cap) + " roomtype " + str(sd[s].roomtype) + " course " +  str(sd[s].course_name) + " enrol=" + str(str(secEnroll[s])) + " room=" + str(str(roomAssign[s,d,t]))
		gtt[s] = table
	return gtt	
		


		
def makeProfTables(flist,x,z,w,u,GDTS):
	[secEnroll, roomAssign, groupsInSec] = createStats(flist,x,z,w,u,GDTS)
	gtt = {}
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	
	list_of_colors = ['\\cellcolor{purple!30}', '\\cellcolor{blue!30}', '\\cellcolor{green!30}', '\\cellcolor{orange!30}', '\\cellcolor{gray!30}', '\\cellcolor{purple!60}', '\\cellcolor{blue!60}', '\\cellcolor{green!60}','\\cellcolor{orange!60}', '\\cellcolor{gray!60}', '\\cellcolor{purple!90}', '\\cellcolor{blue!90}', '\\cellcolor{green!90}','\\cellcolor{orange!90}', '\\cellcolor{gray!90}', '\\cellcolor{pink!30}','\\cellcolor{pink!60}', '\\cellcolor{pink!90}'  ]
	
	
	
	Q = [(p,s) for p in P for s in S if sd[s].course_num in  pd[p].secNums ]
	PS = [(p,s) for p in P for s in S if s in pd[p].secNums]

	for p in P:
		color = {}
		U = [s for s in pd[p].secNums]
		if len(U) > len(list_of_colors):
			print(len(U))
			exit('not enough colors')
			
		count = 0
		for s in U:
			color[s] = list_of_colors[count]
			count+=1
				
				
		table = [['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', '']]
		for d in D:
			for t in T:
				for s in U:
					if sd[s].course_name[2] == '&':
						sd[s].course_name = sd[s].course_name[:2] + sd[s].course_name[3:]
					for r in R:
						if ((s,d,t,r) in z) and (z[s,d,t,r] > 0.5):
							print(s,d,t,r)
							table[d][t] = table[d][t] + color[s] + 'Section ' + str(s) + " Extended "  +str(sd[s].isExtended) + " cap " + str(sd[s].cap) + " roomtype " + str(sd[s].roomtype) + " course " +  str(sd[s].course_name) + " enrol=" + str(str(secEnroll[s])) + " room=" + str(str(roomAssign[s,d,t]))
		gtt[p] = table
	return gtt	
		



def makeGroupTables(flist,x,z,w,u,GDTS):
	[secEnroll, roomAssign, groupsInSec] = createStats(flist,x,z,w,u,GDTS)
	gtt = {}
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	
	list_of_colors = ['\\cellcolor{purple!30}', '\\cellcolor{blue!30}', '\\cellcolor{green!30}', '\\cellcolor{orange!30}', '\\cellcolor{gray!30}', '\\cellcolor{purple!60}', '\\cellcolor{blue!60}', '\\cellcolor{green!60}','\\cellcolor{orange!60}', '\\cellcolor{gray!60}', '\\cellcolor{purple!90}', '\\cellcolor{blue!90}', '\\cellcolor{green!90}','\\cellcolor{orange!90}', '\\cellcolor{gray!90}', '\\cellcolor{pink!30}','\\cellcolor{pink!60}', '\\cellcolor{pink!90}'  ]
	
	
	
	
	for g in G:
		color = {}
		U = [s for s in S if ( ((g,s) in x) and (x[g,s] > 0.5)    )  ]
		if len(U) > len(list_of_colors):
			print(len(U))
			exit('not enough colors')
			
		count = 0
		for s in U:
			color[s] = list_of_colors[count]
			count+=1
				
				
		table = [['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', ''], ['','', '', '', '', '', '', '']]
		for d in D:
			for t in T:
				for s in S:
					if sd[s].course_name[2] == '&':
						sd[s].course_name = sd[s].course_name[:2] + sd[s].course_name[3:]		
					if ((g,d,t,s) in GDTS) and (u[g,d,t,s] > 0.5):
						table[d][t] = table[d][t] + color[s] + 'Sec ' + str(s) + " Ext "  +str(sd[s].isExtended) + " cap " + str(sd[s].cap) + " rt " + str(sd[s].roomtype) + " " +  str(sd[s].course_name) + " enrl=" + str(str(secEnroll[s])) + " rm=" + str(str(roomAssign[s,d,t]))
		gtt[g] = table
	return gtt	



def create_latext_tablesSecs(flist,x,z,w,u,GDTS):
	latex_header = """
		\\documentclass{article} \\usepackage[margin=0.5in]{geometry}
\\usepackage[table]{xcolor} 
\\newcommand{\\textgr}[1]{\\cellcolor{gray!40}\\textbf{#1}}
\\usepackage{array}
\\newcolumntype{z}{>{\\raggedright\\arraybackslash}p{0.75in}}
\\setlength{\\tabcolsep}{0pt}
\\nofiles

\\begin{document}
	"""
	latex_table = """
	
\\begin{tabular}{|c|z|z|z|z|z|z|z|z|}
\\hline
 \\textgr{} & \\textgr{1} & \\textgr{2} & \\textgr{3} & \\textgr{4} & \\textgr{L} & \\textgr{5} & \\textgr{6} & \\textgr{7}\\\\
\\hline
\\textgr{M}  &  %s  &  %s  &  %s &  %s  & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{T}  & %s & %s & %s & %s & %s &  %s & %s & %s \\\\[65pt]
\\hline
\\textgr{W}  &  %s  & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{R} &  %s   & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{F} &  %s & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\end{tabular}

	"""
	latex_footer = """\\end{document}""" 
	[S, sd,P,pd,G,gd,C,cd,roomtypeList, roomMaster, R, rd, roomToType, D,T, Lunch] = plist
	
	
	
#	gtt = {}
	gtt = makeSecTables(plist,x,z,w,u,GDTS)	
	[secEnroll, roomAssign, groupsInSec] = createStats(plist,x,z,w,u,GDTS)
	
	name="sections"
	tf = open(name + '.tex','w')
	tf.write(latex_header)
	
	for s in S:
		table = gtt[s]
		atuple = table_to_tuple(table)
		tf.write(latex_table % atuple)
	#have it print the professers required section list	
		myGroupsEnr = ""
		for g in groupsInSec[s]:
			myGroupsEnr = myGroupsEnr + ", " + str(g) 
			
		groupInfo = "List of groups enrolled is  %s " % myGroupsEnr 
		tf.write('\n' + groupInfo)
		tf.write('\n' + '\\newpage')
		
	tf.write(latex_footer)
	tf.close()
	os.system('pdflatex ' + name + '.tex')







def create_latext_tablesProfs(flist,x,z,w,u,GDTS):
	latex_header = """
		\\documentclass{article} \\usepackage[margin=0.5in]{geometry}
\\usepackage[table]{xcolor} 
\\newcommand{\\textgr}[1]{\\cellcolor{gray!40}\\textbf{#1}}
\\usepackage{array}
\\newcolumntype{z}{>{\\raggedright\\arraybackslash}p{0.75in}}
\\setlength{\\tabcolsep}{0pt}
\\nofiles

\\begin{document}
	"""
	latex_table = """
	
\\begin{tabular}{|c|z|z|z|z|z|z|z|z|}
\\hline
 \\textgr{} & \\textgr{1} & \\textgr{2} & \\textgr{3} & \\textgr{4} & \\textgr{L} & \\textgr{5} & \\textgr{6} & \\textgr{7}\\\\
\\hline
\\textgr{M}  &  %s  &  %s  &  %s &  %s  & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{T}  & %s & %s & %s & %s & %s &  %s & %s & %s \\\\[65pt]
\\hline
\\textgr{W}  &  %s  & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{R} &  %s   & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{F} &  %s & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\end{tabular}

	"""
	latex_footer = """\\end{document}""" 
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	
	
	
#	gtt = {}
	gtt = makeProfTables(flist,x,z,w,u,GDTS)	
	
	
	name="profs"
	tf = open(name + '.tex','w')
	tf.write(latex_header)
	
	for p in P:
		table = gtt[p]
		atuple = table_to_tuple(table)
		tf.write(latex_table % atuple)
	#have it print the professers required section list	
		mycourses = ""
		for c in pd[p].course_numbers_list:
			mycourses = mycourses + ", " + str(c) 
			
		groupInfo = "Num of sections teaching %s, num of periods %s, courses teaching %s " % (str(pd[p].numSections), str(pd[p].numPeriods), mycourses )
		tf.write('\n' + groupInfo)
		tf.write('\n' + '\\newpage')
		
	tf.write(latex_footer)
	tf.close()
	os.system('pdflatex ' + name + '.tex')






def create_latext_tablesGroups(flist,x,z,w,u,GDTS):
	latex_header = """
		\\documentclass{article} \\usepackage[margin=0.5in]{geometry}
\\usepackage[table]{xcolor} 
\\newcommand{\\textgr}[1]{\\cellcolor{gray!40}\\textbf{#1}}
\\usepackage{array}
\\newcolumntype{z}{>{\\raggedright\\arraybackslash}p{0.75in}}
\\setlength{\\tabcolsep}{0pt}
\\nofiles

\\begin{document}
	"""
	latex_table = """
	
\\begin{tabular}{|c|z|z|z|z|z|z|z|z|}
\\hline
 \\textgr{} & \\textgr{1} & \\textgr{2} & \\textgr{3} & \\textgr{4} & \\textgr{L} & \\textgr{5} & \\textgr{6} & \\textgr{7}\\\\
\\hline
\\textgr{M}  &  %s  &  %s  &  %s &  %s  & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{T}  & %s & %s & %s & %s & %s &  %s & %s & %s \\\\[65pt]
\\hline
\\textgr{W}  &  %s  & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{R} &  %s   & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{F} &  %s & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\end{tabular}

	"""
	latex_footer = """\\end{document}""" 
	
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	
	W = set([(g,s) for g in G for s in S if (sd[s].course_num in gd[g].course_num_list) ]  )
	(gr,mini1,mini2,mini3) = CreateMiniGroups(x,G,W,S)
	
	
#	gtt = {}
	gtt = makeGroupTables(flist,x,z,w,u,GDTS)	
	
	name="groups"
	tf = open(name + '.tex','w')
	tf.write(latex_header)
	for grp in range(gr):
		g = mini3[grp][0]
		table = gtt[g]
		atuple = table_to_tuple(table)
		tf.write(latex_table % atuple)
		
		mycourses = ""
		for c in gd[g].course_name_list:
			if c[2] == '&':
				c = c[:2] + c[3:]
			mycourses = mycourses + ", " + str(c) 
			
		groupInfo = "Group "  + gd[g].divname+str(grp) +   ", Size of group  " + str(mini1[grp]) + ", num of periods %s, courses required %s " % (str(gd[g].numPeriods), mycourses )
		tf.write('\n' + groupInfo)
		tf.write('\n' + '\\newpage')
		
	tf.write(latex_footer)
	tf.close()
	os.system('pdflatex ' + name + '.tex')



def create_latext_tablesRooms(flist,x,z,w,u,GDTS):
	latex_header = """
		\\documentclass{article} \\usepackage[margin=0.5in]{geometry}
\\usepackage[table]{xcolor} 
\\newcommand{\\textgr}[1]{\\cellcolor{gray!40}\\textbf{#1}}
\\usepackage{array}
\\newcolumntype{z}{>{\\raggedright\\arraybackslash}p{0.75in}}
\\setlength{\\tabcolsep}{0pt}
\\nofiles

\\begin{document}
	"""
	latex_table = """
	
\\begin{tabular}{|c|z|z|z|z|z|z|z|z|}
\\hline
 \\textgr{} & \\textgr{1} & \\textgr{2} & \\textgr{3} & \\textgr{4} & \\textgr{L} & \\textgr{5} & \\textgr{6} & \\textgr{7}\\\\
\\hline
\\textgr{M}  &  %s  &  %s  &  %s &  %s  & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{T}  & %s & %s & %s & %s & %s &  %s & %s & %s \\\\[65pt]
\\hline
\\textgr{W}  &  %s  & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{R} &  %s   & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\textgr{F} &  %s & %s & %s & %s & %s & %s & %s & %s \\\\[65pt]
\\hline
\\end{tabular}

	"""
	latex_footer = """\\end{document}""" 
	(S,R,RT,C,P,G,L,D,T,Lunch,sd,rd,pd,gd,catalog,divSizeDict,divDict) = flist
	
	
	
#	gtt = {}
	gtt = makeRoomTables(flist,x,z,w,u,GDTS)	
		
	
	name="rooms"
	tf = open(name + '.tex','w')
	tf.write(latex_header)
	
	for r in R:
		table = gtt[r]
		atuple = table_to_tuple(table)
		tf.write(latex_table % atuple)
		
		mycourses = ""
		
			
		groupInfo = "Room Name %s, Room type is %s " % (rd[r].name, rd[r].roomtype)
		tf.write('\n' + groupInfo)
		tf.write('\n' + '\\newpage')
		
	tf.write(latex_footer)
	tf.close()
	os.system('pdflatex ' + name + '.tex')






def csvSaveIns(plist,dname):
	[S, sd,P,pd,G,gd,C,cd,roomtypeList, roomMaster, R, rd, roomToType, D,T, Lunch] = plist
	dname = dname + str(random.randint(1,1000000))
	os.mkdir(dname)
	os.chdir(dname)
	
	f = open('f1roomtypes.csv', mode='w') 
	fcsv = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	fcsv.writerow(['roomtype'])
	for rt in roomtypeList:
		fcsv.writerow([rt])
	f.close()		

#	fcsv.writerow(['John Smith', 'Accounting', 'November'])

	f = open('f2roomnums.csv', mode='w') 
	fcsv = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	fcsv.writerow(['roomnum','roomtype', 'roomcap'])
	for r in R:
		fcsv.writerow([r,rd[r].type, rd[r].cap])
	f.close()	

	

	f = open('f3courses.csv', mode='w') 
	fcsv = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	fcsv.writerow(['coursenum','isextended', 'periods','roomtype'])
	for c in C:
		fcsv.writerow([c,int(cd[c].isExtended), cd[c].periods, cd[c].roomtype])
	f.close()	

	f = open('f4profnums.csv', mode='w') 
	fcsv = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	fcsv.writerow(['profnum'])
	for p in P:
		fcsv.writerow([p])
	f.close()	

#sectionNums, each has a prof, coursenum, cap
	
	f = open('f5sections.csv', mode='w') 
	fcsv = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	fcsv.writerow(['sectionnum','profnum','coursenum','sectioncap'])
	for s in S:
		findProf = [p for p in P if s in pd[p].secs]
		p = findProf[0]
		fcsv.writerow([s,p,sd[s].course_name,sd[s].cap])
	f.close()

#student, numofcourses, list of courses

	f = open('f6students.csv', mode='w') 
	fcsv = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	
	headStudent = ['student','numofcourses', 'c0', 'c1','c2', 'c3', 'c4', 'c5','c6','c7','c8','c9','c10','c11','c12','c13','c14','c15','c16','c17','c18','c19'] 
	fcsv.writerow(headStudent)
	
	for g in G:
		student = [g,gd[g].numCourses] + gd[g].course_name_set
		for i in range(22-len(student)):
			student.append('')
		fcsv.writerow(student)
	f.close()





class Course():
	def __init__(self,name, number, isExtended, periods, roomtype, cap, parent):
		self.number = number
		self.name = name
		self.isExtended = isExtended
		self.periods = periods
		self.roomtype = roomtype
		self.numSections = 0
		self.totalCapCourse = 0
		#what is the capacity of each section, some can be big or small
		self.cap = cap
		self.parent = parent
		self.secList = []
		self.parentCourseName = ''
		self.iAmParent = None
		self.partOfSat = None
		self.partOfFam = None
			
	def create_section(self,minsize, maxsize):
		capOfNewSec = random.randint(minsize,maxsize)
		self.cap[self.numSections] = capOfNewSec
		self.secList.append(self.numSections)	
		self.numSections+=1
		self.totalCapCourse += capOfNewSec
		
		
	def printInfo(self):
		print(self.number,self.isExtended,self.periods, self.roomtype)
		if len(self.cap) > 0:
			for (k,v) in self.cap.items():
				print('sec ', k,'cap ', v, end=' ')	
			print()	


class Section():
	def __init__(self,number,isExtended,periods,roomtype,cap,course_num,course_name,parent,sec_name, labtie, iAmParent):
		self.number = number
		self.load = 0
		self.discp = ''
		self.isExtended = isExtended
		self.periods = periods
		self.roomtype = roomtype
		self.cap = cap
		self.course_num = course_num
		self.course_name = course_name
		self.parent = parent
		self.sec_name = sec_name
		self.labtie = labtie
		self.iAmParent = iAmParent

class Division():
	def __init__(self,name, number,term,divsize,course_name_list,course_num_list):
		self.number = number
		self.size = divsize
		self.name = name
		self.term = term
		self.numCourses = 0
		self.numPeriods = 0
		#a student can not request 2 copies of the same course
		self.course_name_list = course_name_list
		self.course_num_list = course_num_list
		
	def addCourse(self,course):
			self.course_name_set.append(course.number)
			self.numCourses +=1
			self.numPeriods += course.periods


		
	


class Student():
	def __init__(self,number, name, divname, divnumber,course_name_list,course_num_list):
		self.number = number
		self.name = name
		self.divname = divname
		self.divnumber = divnumber
		
		self.numCourses = 0
		self.numPeriods = 0
		#a student can not request 2 copies of the same course
		self.course_name_list = course_name_list
		self.course_num_list = course_num_list
		self.course_num_set = set([])
		
	def addCourse(self,course):
			self.course_name_set.append(course.number)
			self.numCourses +=1
			self.numPeriods += course.periods


				
class Room():
	def __init__(self,name,number,roomtype, cap):
		self.number = number
		self.name = name
		self.roomtype = roomtype
		self.cap = cap
				
				
class Prof():
	def __init__(self,number,discp):
		self.number = number
		self.numSections = 0
		self.courses = []
		self.numPeriods = 0
		self.course_numbers_list = []
		self.secs = []
		self.secNums = []
		self.discp = discp
		self.load = 0
	
	def addSec(self,s):
		self.secs.append(s)
		self.secNums.append(s.number)
		self.numPeriods += s.periods
		self.load+=s.load
		
	
	
	def addSecNum(self,s):
		self.secs.append(s)	
		
	def addCourse(self,course):
		self.courses.append(course)
		self.numSections+=1
		self.numPeriods+=course.periods
		self.course_numbers_list.append(course.number)
		
def makeNewInst(threads,timeblock,curricu):	

	plist = readCSVPandas(curricu)
	(divDict,term2Requests,divSizeDict,catalog,courseNametoNum,rd,roomTypeList,roomTypeNameToNum) = plist
	(gd,sd) = createStudentsAndSections(plist)
	
	plist = (divDict,term2Requests,divSizeDict,catalog,courseNametoNum,rd,roomTypeList,roomTypeNameToNum,gd,sd)
	
	pd = makeProfs(plist)
	
	plist = (divDict,term2Requests,divSizeDict,catalog,courseNametoNum,rd,roomTypeList,roomTypeNameToNum,gd,sd,pd)
	
	flist = prepInstanceForSolvers(plist)
	
	time = timeblock
	

	copyflist = copy.deepcopy(flist)

	(x_m,x2_m) = greedy(copyflist)
	
	(x_m,x2_m,flist) = secCP(threads, flist,time,x_m,x2_m)
	return (x_m,x2_m,flist)
	
	
	
	

def main():	
	
	choice = sys.argv[1].upper()
	if choice == 'TABLES':
		(flist, x_m, z_m, w_m, u_m, GDTS) = pickle.load( open("solX.p", "rb") )
		create_latext_tablesGroups(flist,x_m,z_m,w_m,u_m, GDTS)
		create_latext_tablesProfs(flist,x_m,z_m,w_m,u_m, GDTS)
		create_latext_tablesRooms(flist,x_m,z_m,w_m,u_m, GDTS)
		exit()
	
	
	threads = int(sys.argv[2])
	
	timeblock = int(sys.argv[3])
	
	curricu = sys.argv[4]
	print(curricu)
	
	
	if choice == 'ROUND':
		for k in range(20):
			solD = {}
			#initial runs
			for i in range(1):
				print('Initial ',i,'of ',0)
				tup1 = makeNewInst(threads, timeblock,curricu)
				(x_m,x2_m,flist) = tup1
				f1 = s2(x_m, x2_m, flist,threads,timeblock, None)
				if f1 == None:
					timeblock*=2
					tup1 = makeNewInst(threads, timeblock)
					(x_m, x2_m, flist) = tup1
					f1 = s2(x_m, x2_m, flist, threads, timeblock, None)
				else:
					(h_m,tabu1,status) = f1

				tup1 = (x_m,x2_m,flist)
				solD[i] = copy.deepcopy((tup1,h_m,tabu1,status))



			for i in range(3):
				j = i%1
				print('Continuing ',j,'with i value',i)
				(tup1,h_m,tabu1,status) = solD[j]
				if status < 10:
					print("Problem ",k,"Solved ")
					break
				(x_m,x2_m,flist) = tup1
				(x_m,x2_m,flist) = secCP(threads, flist,timeblock,x_m,x2_m,tabu1)
				f1 = s2(x_m, x2_m, flist,threads,timeblock, h_m)
				if f1 == None:
					timeblock *= 2
					f1 = s2(x_m, x2_m, flist, threads, timeblock, h_m)
					(h_m, tabu1, status) = f1
				else:
					(h_m, tabu1, status) = f1


				tup1 = (x_m,x2_m,flist)
				solD[j] = copy.deepcopy((tup1,h_m,tabu1,status))



	
		
	if choice == 'NEW':
		makeNewInst(threads, timeblock, curricu)  # create new instance and solve the sectioning problem 3 hours
		time.sleep(3)
	
# 	savedSol = 	pickle.load( open("solX.p", "rb") )
	
	if choice == 'EDGE':
		for curricu in ('easy', 'medium', 'medium2', 'hard'):
			print(curricu)
			for i in range(5):
				print('attempt ',i)
				(x_m,x2_m,flist) = makeNewInst(threads, timeblock, curricu )  # create new instance and solve the sectioning problem 3 hours
	# 			time.sleep(3)
	# 			(x_m,x2_m,flist) = pickle.load( open("saveX.p", "rb") )
	#			s2(x_m, x2_m, flist, threads,timeblock)
	
	
	if choice == 'BENCH':
		sectimes = [100,600,1800]
		for j in range(3):
			for curricu in ('easy', 'medium', 'medium2', 'hard'):
				for i in range(3):
					print("attempt ",i,curricu)
					timeblock = sectimes[j]
					(x_m,x2_m,flist) = makeNewInst(threads, timeblock, curricu )  
					timeblock = 600
					s2(x_m, x2_m, flist, threads,timeblock)
	
	
	
	if choice == 'BOTH':
		for i in range(5):
			print("attempt ",i)
			(x_m,x2_m,flist) = makeNewInst(threads, timeblock, curricu )  # create new instance and solve the sectioning problem 3 hours
# 			time.sleep(3)
# 			(x_m,x2_m,flist) = pickle.load( open("saveX.p", "rb") )
			s2(x_m, x2_m, flist, threads,timeblock)

	

	if choice == 'OLD':
		(x_m,x2_m,flist) = pickle.load( open("saveX.p", "rb") )
		s2(x_m, x2_m, flist, threads,timeblock)
	
	
	
	
			

main()	
		
	
