#!/Users/martinmccullagh/anaconda/bin/python
##!/Library/Frameworks/Python.framework/Versions/Current/bin/python
#USAGE : 
import scipy
import sys
import os
import numpy
import math
import MDAnalysis
from MDAnalysis.analysis.align import *

radians_to_degrees = 180.0/3.1415926535
base_atom_select = "name N9 or name N7 or name C8 or name C5 or name C4 or name N3 or name C2 or name N1 or name C6 or name O6 or name N2 or name N6 or name O2 or name N4 or name O4 or name C5M or name N2"
nucleic_res_select = "resname DA or resname DA3 or resname DA5 or resname ADE or resname DT or resname DT3 or resname DT5 or resname THY or resname DG or resname DG3 or resname DG5 or resname GUA or resname DC or resname DC3 or resname DC5 or resname CYT or resname DU or resname DU3 or resname DU5 or resname URA or resname AP3"

# Subroutines
def get_rot_mat(vec1,vec2):
	v = numpy.cross(vec1,vec2)
	cosTheta = numpy.dot(vec1,vec2)
	vx = numpy.matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]],dtype=float)
	rot = numpy.diag([1.0,1.0,1.0]) + vx + numpy.dot(vx,vx)*(1-cosTheta)/math.pow(linalg.norm(v),2)
	return rot

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, traj_file
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='topfile':
				top_file = value
			elif option.lower()=='trajfile':
				traj_file = value
			else :
				print "Option:", option, " is not recognized"

def computePbcDist(r1,r2,box):
	dist = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist += temp*temp

	dist = math.sqrt(dist)
	return dist;

def computeBaseAxes(nucl):


	n_residues = len(nucl.resids())

	axes = numpy.zeros((n_residues,3,3),dtype=float)

	for i in range(0,n_residues):
		# select the residue of interest
		selection = "resid " + str(nucl.resids()[i])
		res_univ = nucl.selectAtoms(selection)
		# get the three atom positions
		n1_pos = res_univ.N1.pos
		n9_pos = res_univ.N9.pos
		n7_pos = res_univ.N7.pos
		# determine the three axes using the three atom positions
		r1 = n1_pos-n9_pos
		r1 /= math.sqrt(numpy.dot(r1,r1))
		t1 = n7_pos-n9_pos
		t1 /= math.sqrt(numpy.dot(t1,t1))
		r3 = numpy.cross(r1,t1)
		r3 /= math.sqrt(numpy.dot(r3,r3))
		r2 = numpy.cross(r3,r1)
		r2 /= math.sqrt(numpy.dot(r2,r2))
		# save the axes for each base
		axes[i][0] = r1
		axes[i][1] = r2
		axes[i][2] = r3

	return axes;

def computeBaseDistances(nucl,base_axes,dist_file_pointer):
	global base_atom_select
	n_residues = len(nucl.resids())
	coms = numpy.zeros((n_residues,3),dtype=float)

	for i in range(0,n_residues):
		# select the residue of interest
		selection = "resid " + str(nucl.resids()[i])
		res_univ = nucl.selectAtoms(selection)
		base_univ = res_univ.selectAtoms(base_atom_select)
		# get the center of mass of the base
		coms[i] = base_univ.centerOfMass()

	for i in range(0,n_residues-1):
		for j in range (i+1,n_residues):
			dist = coms[j]-coms[i]
			slide = numpy.dot(dist,base_axes[i][0])
			shift = numpy.dot(dist,base_axes[i][1])
			rise = numpy.dot(dist,base_axes[i][2])
			dist_file_pointer.write("%3d-%3d %8.3f %8.3f %8.3f\n" % (i+1,j+1,slide,shift,rise))

def computeBaseAngles(nucl,base_axes,ang_file_pointer):
	global radians_to_degrees
	n_residues = len(nucl.resids())

	for i in range(0,n_residues-1):
		for j in range (i+1,n_residues):
			# Roll is the rotation about the major axis (r1) 
			# determine this by first projecting r2' into the r2xr3 plane
			r2p = numpy.dot(base_axes[j][1],base_axes[i][1])*base_axes[i][1] + numpy.dot(base_axes[j][1],base_axes[i][2])*base_axes[i][2]
			r2p /= math.sqrt(numpy.dot(r2p,r2p))
			cos_roll = numpy.dot(r2p,base_axes[i][1])
			sin_roll = numpy.dot(r2p,base_axes[i][2])
			roll = math.atan2(sin_roll,cos_roll)*radians_to_degrees
			# Tilt is the rotation about the in-plane minor axis (r2)
			# we project r3' into the r1xr3 plane
			r3p = numpy.dot(base_axes[j][2],base_axes[i][0])*base_axes[i][0]+numpy.dot(base_axes[j][2],base_axes[i][2])*base_axes[i][2]
			r3p /= math.sqrt(numpy.dot(r3p,r3p))
			cos_tilt = numpy.dot(r3p,base_axes[i][2])
			sin_tilt = numpy.dot(r3p,base_axes[i][0])
			tilt = math.atan2(sin_tilt,cos_tilt)*radians_to_degrees
			# Twist is the rotation about the out-of-plane axis (r3)
			# we project r1' into the r1xr2 plane
			r1p = numpy.dot(base_axes[j][0],base_axes[i][0])*base_axes[i][0]+numpy.dot(base_axes[j][0],base_axes[i][1])*base_axes[i][1]
			r1p /= math.sqrt(numpy.dot(r1p,r1p))
			cos_twist = numpy.dot(r1p,base_axes[i][0])
			sin_twist = numpy.dot(r1p,base_axes[i][1])
			twist = math.atan2(sin_twist,cos_twist)*radians_to_degrees
			# print values to file for this pair of bases
			ang_file_pointer.write("%3d-%3d %8.3f %8.3f %8.3f\n" % (i+1,j+1,roll,tilt,twist))


# Main Program

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Trajectory file:", traj_file

# initiate coordinate universe
coord = MDAnalysis.Universe(top_file, traj_file)

# Select only nucleic atoms
nucl = coord.selectAtoms(nucleic_res_select)

print nucl

# open output files
dist_file_pointer = open("dist.dat",'w')
ang_file_pointer = open("ang.dat",'w')

# Loop through trajectory
for ts in coord.trajectory:

	# compute axes of base
	base_axes = computeBaseAxes(nucl)

	# compute distance parameters for each pair of bases
	computeBaseDistances(nucl,base_axes,dist_file_pointer)

	# compute angle parameters for each pair of bases
	computeBaseAngles(nucl,base_axes,ang_file_pointer)

dist_file_pointer.close
ang_file_pointer.close

# subtract the average from the selected coordinates
#disp_vec = ref_sel.centerOfMass() 
#ref_sel.translate(-disp_vec)
#coord.atoms.translate(-inp_sel.centerOfMass())

#compute rotation matrix
#R, rmsd = MDAnalysis.analysis.align.rotation_matrix(inp_sel.positions,ref_sel.positions)
#print R, disp_vec

#Apply rotation and translation and then print new pdb
#coord.atoms.rotate(R)
#coord.atoms.translate(disp_vec)
#coord.atoms.write(out_pdb_file)
