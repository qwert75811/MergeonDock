# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 16:51:19 2024

@author: Xhamrock Studio
"""

from pymol import cmd, cgo

class PyMOLGridBox:
    def __init__(self, center=[0, 0, 0], size=[20, 20, 20], space=0.375, alpha=0.5):
        
        self.center = center
        self.size = size
        self.space = space
        self.alpha = alpha  # 添加 alpha 透明度变量
        self.colors = [
            [0, 0, 1],  # Front face color: red
            [0, 0, 1],  # Back face color: green
            [1, 0, 0],  # Left face color: blue
            [1, 0, 0],  # Right face color: yellow
            [0, 1, 0],  # Top face color: cyan
            [0, 1, 0]   # Bottom face color: magenta
        ]
        cmd.set('auto_zoom', 'off')
        self.draw_colored_box()
        

    def draw_colored_box(self):
        x, y, z = self.center
        sizex, sizey, sizez = self.size
        half_sizex = (sizex * self.space) / 2.0
        half_sizey = (sizey * self.space) / 2.0
        half_sizez = (sizez * self.space) / 2.0

        box = []

        # Front face
        box.extend([cgo.COLOR, *self.colors[0], cgo.ALPHA, self.alpha,
                    cgo.BEGIN, cgo.TRIANGLES,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z + half_sizez,
                    cgo.END])

        # Back face
        box.extend([cgo.COLOR, *self.colors[1], cgo.ALPHA, self.alpha,
                    cgo.BEGIN, cgo.TRIANGLES,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z - half_sizez,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z - half_sizez,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z - half_sizez,
                    cgo.END])

        # Left face
        box.extend([cgo.COLOR, *self.colors[2], cgo.ALPHA, self.alpha,
                    cgo.BEGIN, cgo.TRIANGLES,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z - half_sizez,
                    cgo.END])

        # Right face
        box.extend([cgo.COLOR, *self.colors[3], cgo.ALPHA, self.alpha,
                    cgo.BEGIN, cgo.TRIANGLES,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z - half_sizez,
                    cgo.END])

        # Top face
        box.extend([cgo.COLOR, *self.colors[4], cgo.ALPHA, self.alpha,
                    cgo.BEGIN, cgo.TRIANGLES,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z - half_sizez,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x - half_sizex, y + half_sizey, z - half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y + half_sizey, z - half_sizez,
                    cgo.END])

        # Bottom face
        box.extend([cgo.COLOR, *self.colors[5], cgo.ALPHA, self.alpha,
                    cgo.BEGIN, cgo.TRIANGLES,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x - half_sizex, y - half_sizey, z - half_sizez,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z + half_sizez,
                    cgo.VERTEX, x + half_sizex, y - half_sizey, z - half_sizez,
                    cgo.END])
        
        
        cmd.delete('gridbox')
        return box

    def update_center(self, x, y, z):
        self.center = [x, y, z]
        self.draw_colored_box()

    def update_size(self, sizex, sizey, sizez):
        self.size = [sizex, sizey, sizez]
        self.draw_colored_box()

    def update_space(self, space):
        self.space = space
        self.draw_colored_box()





   
