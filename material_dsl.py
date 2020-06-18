"""
Create Blender materials from Python:
    from material_dsl import (BsdfPrincipled, OutputMaterial, Value, VertexColor, make_material)
    color_input = VertexColor(location=(-400, 300), layer_name="my colormap")
    shader = BsdfPrincipled(location=(0, 300), base_color=color_input.color)
    output_material = OutputMaterial(location=(400, 300), surface=shader.BSDF)
    
    make_material("colored mesh", output_material)
    
The set of shader nodes needs to be expanded, but this looks a lot better than
adding nodes and links by hand.
"""

from __future__ import annotations
import bpy

import functools
from dataclasses import dataclass
from typing import List, Any, Dict, Union, Tuple


def decorator(f):
    """Creates a paramatric decorator from a function. The resulting decorator
    will optionally take keyword arguments."""
    @functools.wraps(f)
    def decoratored_function(*args, **kwargs):
        if args and len(args) == 1:
            return f(*args, **kwargs)

        if args:
            raise TypeError(
                "This decorator only accepts extra keyword arguments.")

        return lambda g: f(g, **kwargs)

    return decoratored_function


@dataclass
class Promise:
    graph: Graph
    output: Output


@dataclass
class Graph:
    nodes: List[Node]
    links: List[Tuple[Output, Input]]

    @property
    def root(self):
        return self.nodes[0]

    def __getattr__(self, name):
        return Promise(self, Output(self.root, name))


@dataclass
class Output:
    node: Node
    name: str


@dataclass
class Input:
    node: Node
    name: str


@dataclass
class Value:
    value: Any


@dataclass
class Node:
    name: str
    properties: Dict[str, Any]
    input_defaults: Dict[Union[int, str], Union[Value, Output]]


@decorator
def node(f, properties=["location"]):
    @functools.wraps(f)
    def g(*args, **kwargs):
        name = f.__name__
        property_values = {}
        input_defaults = {}

        links = []
        nodes = [Node(name, property_values, input_defaults)]

        def merge_graph(g):
            for n in g.nodes:
                if n not in nodes:
                    nodes.append(n)
            for link in g.links:
                if link not in links:
                    links.append(link)

        for i, a in enumerate(args):
            if isinstance(a, Value):
                input_defaults[i] = a
            elif isinstance(a, Promise):
                merge_graph(a.graph)
                links.append((a.output, Input(nodes[0], i)))

        for k, v in kwargs.items():
            if k in properties:
                property_values[k] = v
            elif isinstance(v, Value):
                input_defaults[k] = v
            elif isinstance(v, Promise):
                merge_graph(v.graph)
                links.append((v.output, Input(nodes[0], k)))

        return Graph(nodes, links)
    return g


@node(properties=["location"])
def BsdfPrincipled(**kwargs):
    pass


@node(properties=["location"])
def OutputMaterial(**kwargs):
    pass


@node(properties=["location", "layer_name"])
def VertexColor(**kwargs):
    pass


@node(properties=["location"])
def MixShader(*args, **kwargs):
    pass


@node
def BsdfTransparent(**kwargs):
    pass


def demangle(name):
    if isinstance(name, int):
        return name

    def cap(s):
        return s[0].upper() + s[1:]

    return ' '.join([cap(w) for w in name.split('_')])


def make_material(name: str, graph: Graph):
    material = bpy.data.materials.new(name)
    material.use_nodes = True
    nodes = material.node_tree.nodes
    nodes.clear()

    nodemap = {}
    for n in graph.nodes:
        s = nodes.new(type=f"ShaderNode{n.name}")
        nodemap[id(n)] = s
        for k, v in n.properties.items():
            setattr(s, k, v)
        for k, v in n.input_defaults.items():
            key = demangle(k)
            s.inputs[key].default_value = v.value

    links = material.node_tree.links
    for (o, i) in graph.links:
        node_out = nodemap[id(o.node)]
        node_in = nodemap[id(i.node)]
        links.new(node_out.outputs[demangle(o.name)],
                  node_in.inputs[demangle(i.name)])

    return material