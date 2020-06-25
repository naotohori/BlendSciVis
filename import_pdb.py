from Bio.PDB import PDBParser
import bpy
import bmesh
#from mathutils import Matrix, Vector
import mathutils

objects = bpy.data.objects
meshes = bpy.data.meshes
collection = bpy.data.collections[0]

# Delete existing objects
for obj in bpy.data.objects:
    #if obj.name[0:11] == 'atom_sphere':
    #    objects.remove(objects[obj.name])
    if obj.name not in ('Camera', 'Light'):
        objects.remove(objects[obj.name])

def put_sphere(col, name, coord, rad):
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

    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.ops.object.shade_smooth()

parser = PDBParser()
structure = parser.get_structure('azo1000', '1000.pdb')

def clear_collection(name):
    coll = bpy.data.collections.get(name)
    if coll:
        #remove_collection_objects = True
        if True:
            obs = [o for o in coll.objects if o.users == 1]
            while obs:
                bpy.data.objects.remove(obs.pop())

        bpy.data.collections.remove(coll)

# Delete exisint collections except default "Collection"
for col in bpy.data.collections:
    if col.name != 'Collection':
        clear_collection(col.name)

# Generate collections
RNA = bpy.data.collections.new("RNA")
Mg = bpy.data.collections.new("Mg")
K = bpy.data.collections.new("K")
Cl = bpy.data.collections.new("Cl")
bpy.context.scene.collection.children.link(RNA)
bpy.context.scene.collection.children.link(Mg)
bpy.context.scene.collection.children.link(K)
bpy.context.scene.collection.children.link(Cl)

atom_id = 0
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                atom_id += 1    
                if atom_id > 196:
                    continue
                #if atom_id < 10:
                #    print('model: %i' % model.get_id())
                #    print('chain: %s' % chain.get_id())
                #    print('residue: %i' % residue.get_id()[1])
                #    print('atom: %s' % atom.get_id())
                #    # access to atom
                #    # structure[0]['R'][196]['P']
                #atom = structure[0]['A'][1]['CA']
                name = 'sphere_%03i' % atom_id
                put_sphere(RNA, name, atom.get_coord(), 1.0)
