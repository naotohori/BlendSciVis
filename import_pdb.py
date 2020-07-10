
FLG_DEBUG = False
DEBUG_NUM_ATOM_PER_TYPE = 20

FLG_SUBSURF = False
SUBSURF_LEVEL = 1
FLG_SMOOTH = False

from Bio.PDB import PDBParser
import math
import bpy
import bmesh
from mathutils import Matrix, Vector, Euler

#import sys
#import os
## https://blender.stackexchange.com/questions/51044/how-to-import-a-blender-python-script-in-another
#dir = os.path.dirname(bpy.data.filepath)
#if not dir in sys.path:
#    sys.path.append(dir )
#
#import put_objects
#
## this next part forces a reload in case you edit the source after you first start the blender session
#import imp
#imp.reload(put_objects)
#
## this is optional and allows you to call the functions without specifying the package name
#from put_objects import *

objects = bpy.data.objects
meshes = bpy.data.meshes
collection = bpy.data.collections[0]
materials = bpy.data.materials

# Delete existing objects
for obj in objects:
    #if obj.name[0:11] == 'atom_sphere':
    #    objects.remove(objects[obj.name])
    if obj.name[0:5] == 'Light':
        continue
    if obj.name in ('Camera',):
        continue

    objects.remove(objects[obj.name])

def put_sphere(col, name, coord, rad, mat):
    '''
    Adopted from 
    https://blender.stackexchange.com/questions/93298/create-a-uv-sphere-object-in-blender-from-python
    '''

    # Create an empty mesh and the object.
    mesh = meshes.new(name)
    basic_sphere = objects.new(name, mesh)

    # Add the object into the scene.
    col.objects.link(basic_sphere)

    # Select the newly created object
    bpy.context.view_layer.objects.active = basic_sphere
    basic_sphere.select_set(True)

    # Construct the bmesh sphere and assign it to the blender mesh.
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(bm, u_segments=16, v_segments=8, diameter=2*rad)
    bm.to_mesh(mesh)
    bm.free()

    basic_sphere.location = coord

    if mat is not None:
        if basic_sphere.data.materials:
            basic_sphere.data.materials[0] = mat
        else:
            basic_sphere.data.materials.append(mat)

    if FLG_SUBSURF:
        bpy.ops.object.modifier_add(type='SUBSURF')
        bpy.ops.object.subdivision_set(level=SUBSURF_LEVEL)
    if FLG_SMOOTH:
        bpy.ops.object.shade_smooth()

parser = PDBParser()
structure = parser.get_structure('azo1000', '1000.pdb')

# Delete exisint collections except default "Collection"
def clear_collection(name):
    coll = bpy.data.collections.get(name)
    if coll:
        #remove_collection_objects = True
        if True:
            obs = [o for o in coll.objects if o.users == 1]
            while obs:
                bpy.data.objects.remove(obs.pop())

        bpy.data.collections.remove(coll)

for col in bpy.data.collections:
    if col.name != 'Collection':
        clear_collection(col.name)

# Generate collections
RNA = bpy.data.collections.new("RNA")
Mg = bpy.data.collections.new("Mg")
K = bpy.data.collections.new("K")
Cl = bpy.data.collections.new("Cl")
RNAbonds = bpy.data.collections.new("RNAbonds")
bpy.context.scene.collection.children.link(RNA)
bpy.context.scene.collection.children.link(RNAbonds)
bpy.context.scene.collection.children.link(Mg)
bpy.context.scene.collection.children.link(K)
bpy.context.scene.collection.children.link(Cl)


# Delete existing matrials
for mat in materials:
    materials.remove(materials[mat.name])

# Make material
def make_material(name, color_tuple):
    M = materials.new(name)
    
    M.use_nodes = True
    nodes = M.node_tree.nodes
    links = M.node_tree.links
    nodes.clear()  # Clear default nodes
    
    # Create an output for the shader
    node_output = nodes.new(type='ShaderNodeOutputMaterial')
    node_output.location = 400, 300
    
    shader = nodes.new(type='ShaderNodeBsdfPrincipled')
    shader.location = 0, 300 # Location in the node window
    shader.inputs['Base Color'].default_value = color_tuple
    links.new(shader.outputs['BSDF'], node_output.inputs['Surface'])

    return M

def make_emission(name, color_tuple):
    M = materials.new(name)
    
    M.use_nodes = True
    nodes = M.node_tree.nodes
    links = M.node_tree.links
    nodes.clear()  # Clear default nodes
    
    # Create an output for the shader
    node_output = nodes.new(type='ShaderNodeOutputMaterial')
    node_output.location = 400, 300
    
    #shader = nodes.new(type='ShaderNodeBsdfPrincipled')
    #shader.location = 0, 300 # Location in the node window
    #shader.inputs['Base Color'].default_value = color_tuple
    #links.new(shader.outputs['BSDF'], node_output.inputs['Surface'])
    shader = nodes.new(type='ShaderNodeEmission')
    shader.location = 0, 300 # Location in the node window
    shader.inputs['Color'].default_value = color_tuple
    shader.inputs['Strength'].default_value = 75
    links.new(shader.outputs['Emission'], node_output.inputs['Surface'])

    return M

def mp_to_nt(imp):
    return imp // 3 + 1

def azo_segment(imp):
    nt = mp_to_nt(imp)
    if nt < 26:
        return 'P2'
    elif nt in range(30,36) or nt in range(127,133):
        return 'P3'
    elif nt in range(77,83) or nt in range(40,46):
        return 'P4'
    elif nt in range(49, 74):
        return 'P5'
    elif nt in range(83, 114):
        return 'P6'
    elif nt in range(118, 124) or nt in range(162, 168):
        return 'P7'
    elif nt in range(133, 156):
        return 'P8'
    elif nt in range(168, 197):
        return 'P9'
    elif nt in range(156, 162):
        return 'J87'
    elif nt in range(114, 118):
        return 'TH1'
    elif nt in range(36, 40):
        return 'TH2'
    return 'other'


## RNA
RNAmat   = make_material('RNA_mat', (0.35,0.35,0.35,1))
RNAP2mat = make_material('RNAP2_mat', (1,0.5,0,1))
RNAP3mat = make_material('RNAP3_mat', (0.02,0.38,0.57,1))
RNAP4mat = make_material('RNAP4_mat', (0,0.9,0.5,1))
RNAP5mat = make_material('RNAP5_mat', (0,1,0,1))
RNAP6mat = make_material('RNAP6_mat', (0.5,0.9,0.4,1))
RNAP7mat = make_material('RNAP7_mat', (0.25,0.75,0.75,1))
RNAP8mat = make_material('RNAP8_mat', (0,0,1,1))
RNAP9mat = make_material('RNAP9_mat', (0.45,0,0.9,1))
RNAJ87mat = make_material('RNAJ87_mat', (1,0,0,1))
RNATH1mat = make_material('RNATH1_mat', (1,0.6,0.6,1))
RNATH2mat = make_material('RNATH2_mat', (0.98,0,0,1))

segment_to_mat = {'other': RNAmat,
                    'P2': RNAP2mat,
                    'P3': RNAP3mat,
                    'P4': RNAP4mat,
                    'P5': RNAP5mat,
                    'P6': RNAP6mat,
                    'P7': RNAP7mat,
                    'P8': RNAP8mat,
                    'P9': RNAP9mat,
                    'J87': RNAJ87mat,
                    'TH1': RNATH1mat,
                    'TH2': RNATH2mat,
                 }

## RNA bond
RNAbondmat = make_material('RNAbond_mat', (0.8,0.8,0.8,1))

## Ions
MgEmissionmat = make_emission('Emission_mat', (1, 0.15, 0., 1))
Mgmat = make_material('Mg_mat', (1, 0.15, 0., 1))
Kmat = make_material('K_mat', (0.08, 0.26, 0.07, 1))
Clmat = make_material('Cl_mat', (0.12, 0.23, 0.21, 1))

#mesh.materials.append( mat )

#print('model: %i' % model.get_id())
#print('chain: %s' % chain.get_id())
#print('residue: %i' % residue.get_id()[1])
#print('atom: %s' % atom.get_id())
## access to atom
## structure[0]['R'][196]['P']
#atom = structure[0]['A'][1]['CA']

print("Constructing RNA")
atom_id = 0
for residue in structure[0]['R']:
    for atom in residue:
        atom_id += 1
        if FLG_DEBUG:
            if atom_id > DEBUG_NUM_ATOM_PER_TYPE:
                continue

        segment = azo_segment(atom_id)
        print ('%i %s' % (mp_to_nt(atom_id), segment))
        name = 'R%03i' % atom_id
        put_sphere(RNA, name, atom.get_coord(), 0.75, segment_to_mat[segment])

print("Constructing Mg")
atom_id = 0
for residue in structure[0]['M']:
    for atom in residue:
        atom_id += 1
        if FLG_DEBUG:
            if atom_id > DEBUG_NUM_ATOM_PER_TYPE:
                continue
        name = 'Mg%03i' % atom_id

        #if atom_id == 28:
        if atom_id == 76:
            put_sphere(Mg, name, atom.get_coord(), 0.6, MgEmissionmat)
        else:
            put_sphere(Mg, name, atom.get_coord(), 0.6, Mgmat)

print("Constructing Cl")
atom_id = 0
for residue in structure[0]['C']:
    for atom in residue:
        atom_id += 1
        if FLG_DEBUG:
            if atom_id > DEBUG_NUM_ATOM_PER_TYPE:
                continue
        name = 'Cl%03i' % atom_id
        put_sphere(Cl, name, atom.get_coord(), 0.4, Kmat)

print("Constructing K")
atom_id = 0
for residue in structure[0]['K']:
    for atom in residue:
        atom_id += 1
        if FLG_DEBUG:
            if atom_id > DEBUG_NUM_ATOM_PER_TYPE:
                continue
        name = 'K%03i' % atom_id
        put_sphere(K, name, atom.get_coord(), 0.4, Clmat)

### Bond ####
def put_cylinder(col, name, coord_i, coord_j, rad, mat=None):

    # Create an empty mesh and the object.
    mesh = meshes.new(name)
    cylinder = objects.new(name, mesh)

    # Add the object into the scene.
    col.objects.link(cylinder)

    # Select the newly created object
    bpy.context.view_layer.objects.active = cylinder
    cylinder.select_set(True)

    #  Calculate translation and rotation
    dv = coord_i - coord_j
    length = dv.length
    centre = 0.5 * (coord_i + coord_j)
    phi = math.atan2(dv.y, dv.x) 
    theta = math.acos(dv.z/length) 

    # Construct the bmesh sphere and assign it to the blender mesh.
    bm = bmesh.new()
    #bmesh.ops.create_uvsphere(bm, u_segments=16, v_segments=8, diameter=2*rad)
    bmesh.ops.create_cone(bm, cap_ends=False, 
                                cap_tris=False, segments=12, diameter1=rad, diameter2=rad,
                                 depth=length, #matrix = mtx, calc_uvs=False
                                 )
    bm.to_mesh(mesh)
    bm.free()

    cylinder.location = centre
    bpy.context.object.rotation_euler[1] = theta 
    bpy.context.object.rotation_euler[2] = phi 

    #bpy.ops.object.modifier_add(type='SUBSURF')
    #bpy.ops.object.shade_smooth()
    if mat is not None:
        if cylinder.data.materials:
            cylinder.data.materials[0] = mat
        else:
            cylinder.data.materials.append(mat)

    if FLG_SUBSURF:
        bpy.ops.object.modifier_add(type='SUBSURF')
        bpy.ops.object.subdivision_set(level=SUBSURF_LEVEL)
    if FLG_SMOOTH:
        bpy.ops.object.shade_smooth()


print("Adding bonds")

flg_bond = False
bonds = []
for l in open('azo.psf'):
    if l.find('NBOND') != -1:
        flg_bond = True
        continue
    if flg_bond:
        lsp = l.split()
        for i in range(0,len(lsp),2):
            bonds.append((int(lsp[i]), int(lsp[i+1])))

for b in bonds:
    i = b[0]
    j = b[1]
    name = 'bond_%i_%i' % (i,j)

    if FLG_DEBUG:
        if i > DEBUG_NUM_ATOM_PER_TYPE or j > DEBUG_NUM_ATOM_PER_TYPE:
            continue
    
    name_i = 'R%03i' % i
    name_j = 'R%03i' % j
    if name_i not in RNA.objects:
        continue
    if name_j not in RNA.objects:
        continue
    coord_i = RNA.objects[name_i].location
    coord_j = RNA.objects[name_j].location
        
    #print(coord_i)
    #print(coord_j)
    #print((coord_i - coord_j).length)
    
    put_cylinder(RNAbonds, name, coord_i, coord_j, rad=0.5, mat=RNAbondmat)

print("Render")

scene = bpy.context.scene
scene.render.image_settings.file_format = 'PNG'

for emission in (100, 50, 20, 5):

    materials['Emission_mat'].node_tree.nodes['Emission'].inputs['Strength'].default_value = emission
    
    #bpy.ops.render.render(use_viewport=True)
    #bpy.ops.render.view_show()
    #bpy.ops.image.save_as('Emission%03i.png' % emission)

    scene.render.filepath = "./Emission%03i.png" % emission
    bpy.ops.render.render(write_still = 1)