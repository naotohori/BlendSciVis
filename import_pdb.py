from Bio.PDB import PDBParser
import math
import bpy
import bmesh
from mathutils import Matrix, Vector, Euler

FLG_DEBUG = False

objects = bpy.data.objects
meshes = bpy.data.meshes
collection = bpy.data.collections[0]

# Delete existing objects
for obj in bpy.data.objects:
    #if obj.name[0:11] == 'atom_sphere':
    #    objects.remove(objects[obj.name])
    if obj.name not in ('Camera', 'Light'):
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


    #bpy.ops.object.modifier_add(type='SUBSURF')
    #bpy.ops.object.shade_smooth()

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

# Make material
RNA_material_name = "RNA_mat"
RNAmat = bpy.data.materials.new(RNA_material_name)

RNAmat.use_nodes = True
nodes = RNAmat.node_tree.nodes
links = RNAmat.node_tree.links

# Clear default nodes
nodes.clear()

# Create an output for the shader
node_output = nodes.new(type='ShaderNodeOutputMaterial')
node_output.location = 400, 300

shader = nodes.new(type='ShaderNodeBsdfPrincipled')
shader.location = 0, 300 # Location in the node window
shader.inputs['Base Color'].default_value = (1,0,0,1)
links.new(shader.outputs['BSDF'], node_output.inputs['Surface'])

#mesh.materials.append( mat )

#print('model: %i' % model.get_id())
#print('chain: %s' % chain.get_id())
#print('residue: %i' % residue.get_id()[1])
#print('atom: %s' % atom.get_id())
## access to atom
## structure[0]['R'][196]['P']
#atom = structure[0]['A'][1]['CA']

atom_id = 0
for residue in structure[0]['R']:
    for atom in residue:
        atom_id += 1
        if FLG_DEBUG:
            if atom_id > 10:
                continue

        name = 'R%03i' % atom_id
        put_sphere(RNA, name, atom.get_coord(), 0.5, RNAmat)

atom_id = 0
for residue in structure[0]['M']:
    for atom in residue:
        atom_id += 1
        name = 'Mg%03i' % atom_id
        put_sphere(Mg, name, atom.get_coord(), 0.6, None)

atom_id = 0
for residue in structure[0]['C']:
    for atom in residue:
        atom_id += 1
        name = 'Cl%03i' % atom_id
        put_sphere(Cl, name, atom.get_coord(), 0.4, None)

atom_id = 0
for residue in structure[0]['K']:
    for atom in residue:
        atom_id += 1
        name = 'K%03i' % atom_id
        put_sphere(K, name, atom.get_coord(), 0.4, None)

### Bond ####
def put_cylinder(col, name, coord_i, coord_j, rad):

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
    
    put_cylinder(RNAbonds, name, coord_i, coord_j, rad=0.5)