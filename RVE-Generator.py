# RVE-Generator Porenschluss
# Institut fuer Bildsame Formgebung, RWTH Aachen
# Autor: Xinyang Li, E-Mail: xinyang.li@rwth-aachen.de
# Version 0.1, 13.04.2018

#-------------------------------------------------------------------------------
#Initialisierung
#-------------------------------------------------------------------------------

#Packages importieren
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqusConstants import *
import regionToolset
import os
import math
import re
model = mdb.models[mdb.models.keys()[0]]

#Class definieren
class RVE:
    name = ''
    dimension = '3D' #2D, 3D
    laenge_x = 0.0
    laenge_y = 0.0
    laenge_z = 0.0
    typ_Pore = 'Ellipsoid' #Ellipsoid, Quader, Zylinder
    porenparameter_x = 0.0
    porenparameter_y = 0.0
    porenparameter_z = 0.0
    porenparameter_rx = 0.0
    porenparameter_ry = 0.0
    porenparameter_rz = 0.0
    def __init__(self,name,dimension,typ_Pore,**para):
        self.name = name
        self.dimension = dimension
        self.laenge_x = para['laenge_x']
        self.laenge_y = para['laenge_y']
        if (self.dimension == '3D'):
            self.laenge_z = para['laenge_z']
        self.typ_Pore = typ_Pore
        self.porenparameter_x = para['porenparameter_x']
        self.porenparameter_y = para['porenparameter_y']
        if (self.dimension == '3D'):
            self.porenparameter_z = para['porenparameter_z']
            self.porenparameter_rx = para['porenparameter_rx']
            self.porenparameter_ry = para['porenparameter_ry']
        self.porenparameter_rz = para['porenparameter_rz']
        #Error melden
        if (self.typ_Pore == 'Ellipsoid'):
            if (self.porenparameter_rx != 0.0 or self.porenparameter_ry != 0.0 or self.porenparameter_rz != 0.0):
                raise RuntimeError('Infolge der symmetrischen Einstellung des Modells darf Ellipsoid nicht rotiert werden!')
        elif (self.typ_Pore == 'Quader'):
            if (self.porenparameter_rx != 0.0 and self.porenparameter_y != self.porenparameter_z):
                raise RuntimeError('Die Geometrie- und Rotationsdaten des Quaders passen der symmetrischen Einstellung nicht!')
            if (self.porenparameter_ry != 0.0 and self.porenparameter_x != self.porenparameter_z):
                raise RuntimeError('Die Geometrie- und Rotationsdaten des Quaders passen der symmetrischen Einstellung nicht!')
            if (self.porenparameter_rz != 0.0 and self.porenparameter_x != self.porenparameter_y):
                raise RuntimeError('Die Geometrie- und Rotationsdaten des Quaders passen der symmetrischen Einstellung nicht!')
            if (self.porenparameter_rx % 45 != 0 or self.porenparameter_ry % 45 != 0 or self.porenparameter_rz % 45 != 0):
                raise RuntimeError('Die Rotationswinkeln des Quaders passen der symmetrischen Einstellung nicht!')
        elif (self.typ_Pore == 'Zylinder'):
            if (self.porenparameter_rx % 90 != 0 or self.porenparameter_ry % 90 != 0 or self.porenparameter_rz % 90 != 0):
                raise RuntimeError('Die Rotationswinkeln des Zylinders passen der symmetrischen Einstellung nicht!')
    def sketch_und_part(self):
        if (self.dimension == '3D'):
            #Sketch Wuerfel zeichnen
            self.sketch_Wuerfel = model.ConstrainedSketch(
                name='Seitenansicht_Wuerfel',
                sheetSize=200.0)
            self.sketch_Wuerfel.rectangle(
                point1=(0.0, 0.0),
                point2=(self.laenge_x/2.0, self.laenge_y/2.0)) #x- und y-Symmetrie
            #Part Wuerfel generieren
            self.part_Wuerfel = model.Part(
                name=self.name+'_Wuerfel',
                dimensionality=THREE_D,
                type=DEFORMABLE_BODY)
            self.part_Wuerfel.BaseSolidExtrude(
                sketch=self.sketch_Wuerfel,
                depth=self.laenge_z/2.0) #z-Symmetrie
            #Sketch Pore zeichnen
            self.sketch_Pore = model.ConstrainedSketch(
                name='Seitenansicht_Pore',
                sheetSize=200.0)
            if (self.typ_Pore == 'Ellipsoid' ):
                if (self.porenparameter_x == self.porenparameter_z):
                    self.sketch_Pore.ConstructionLine(
                        point1=(0.0, -100.0),
                        point2=(0.0, 100.0))
                    self.sketch_Pore.EllipseByCenterPerimeter(
                        center=(0.0, 0.0),
                        axisPoint1=(self.porenparameter_x/2.0, 0.0),
                        axisPoint2=(0.0, self.porenparameter_y/2.0))
                    self.sketch_Pore.autoTrimCurve(
                        curve1=self.sketch_Pore.geometry[3],
                        point1=(-self.porenparameter_x/2.0, 0.0))
                    self.sketch_Pore.Line(
                        point1=(0.0, self.porenparameter_y/2.0),
                        point2=(0.0, -self.porenparameter_y/2.0))
                else:
                    #Hilfsgeometrie (3 Oberflaechen), auf den die Ellipse erstellt werden
                    self.sketch_Pore_Skelett = model.ConstrainedSketch(
                        name='Seitenansicht_Pore_Skelett',
                        sheetSize=200.0)
                    self.sketch_Pore_Skelett.rectangle(
                        point1=(0.0, 0.0),
                        point2=(self.porenparameter_x/2.0*3.0, self.porenparameter_y/2.0*3.0))
                    self.part_Pore_Skelett = model.Part(
                        name=self.name+'_Pore_Skelett',
                        dimensionality=THREE_D,
                        type=DEFORMABLE_BODY)
                    self.part_Pore_Skelett.BaseSolidExtrude(
                        sketch=self.sketch_Pore_Skelett,
                        depth=self.porenparameter_z/2.0*3.0)
                    self.transform = self.part_Pore_Skelett.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett.faces[5],
                        sketchUpEdge=self.part_Pore_Skelett.edges[2],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_unvollkommen = model.ConstrainedSketch(
                        name='Seitenansicht_Pore_unvollkommen',
                        sheetSize = 200.0,
                        transform = self.transform)
                    del self.transform
                    self.sketch_Pore_unvollkommen.rectangle(
                        point1=(0.0, 0.0),
                        point2=(-self.porenparameter_x/2.0*2.0, self.porenparameter_y/2.0*2.0))
                    self.part_Pore_Skelett.CutExtrude(
                        sketchPlane=self.part_Pore_Skelett.faces[5],
                        sketchUpEdge=self.part_Pore_Skelett.edges[2],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_unvollkommen,
                        depth=self.porenparameter_z/2.0*2.0,
                        flipExtrudeDirection=OFF)
                    del self.sketch_Pore_unvollkommen
                    #Ellipse XY
                    self.transform = self.part_Pore_Skelett.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett.faces[8],
                        sketchUpEdge=self.part_Pore_Skelett.edges[12],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_Ellipse_XY = model.ConstrainedSketch(
                        name='Seitenansicht_Pore_Ellipse_XY',
                        sheetSize = 200.0,
                        transform = self.transform)
                    del self.transform
                    self.sketch_Pore_Ellipse_XY.EllipseByCenterPerimeter(
                        center=(0.0, 0.0),
                        axisPoint1=(-self.porenparameter_x/2.0, 0.0),
                        axisPoint2=(0.0, self.porenparameter_y/2.0))
                    self.part_Pore_Skelett.Wire(
                        sketchPlane=self.part_Pore_Skelett.faces[8],
                        sketchUpEdge=self.part_Pore_Skelett.edges[12],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_Ellipse_XY)
                    del self.sketch_Pore_Ellipse_XY
                    #Ellipse XZ
                    self.transform = self.part_Pore_Skelett.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett.faces[4],
                        sketchUpEdge=self.part_Pore_Skelett.edges[10],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_Ellipse_XZ = model.ConstrainedSketch(
                        name='Seitenansicht_Pore_Ellipse_XZ',
                        sheetSize = 200.0,
                        transform = self.transform)
                    del self.transform
                    self.sketch_Pore_Ellipse_XZ.EllipseByCenterPerimeter(
                        center=(0.0, 0.0),
                        axisPoint1=(-self.porenparameter_x/2.0, 0.0),
                        axisPoint2=(0.0, self.porenparameter_z/2.0))
                    self.part_Pore_Skelett.Wire(
                        sketchPlane=self.part_Pore_Skelett.faces[4],
                        sketchUpEdge=self.part_Pore_Skelett.edges[10],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_Ellipse_XZ)
                    del self.sketch_Pore_Ellipse_XZ
                    #Ellipse YZ
                    self.transform = self.part_Pore_Skelett.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett.faces[3],
                        sketchUpEdge=self.part_Pore_Skelett.edges[13],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_Ellipse_YZ = model.ConstrainedSketch(
                        name='Seitenansicht_Pore_Ellipse_YZ',
                        sheetSize = 200.0,
                        transform = self.transform)
                    del self.transform
                    self.sketch_Pore_Ellipse_YZ.EllipseByCenterPerimeter(
                        center=(0.0, 0.0),
                        axisPoint1=(-self.porenparameter_y/2.0, 0.0),
                        axisPoint2=(0.0, self.porenparameter_z/2.0))
                    self.part_Pore_Skelett.Wire(
                        sketchPlane=self.part_Pore_Skelett.faces[3],
                        sketchUpEdge=self.part_Pore_Skelett.edges[13],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_Ellipse_YZ)
                    del self.sketch_Pore_Ellipse_YZ
                    #Hilfsgeometrie loeschen (XY-Ebene)
                    self.transform = self.part_Pore_Skelett.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett.faces[8],
                        sketchUpEdge=self.part_Pore_Skelett.edges[24],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_unvollkommen = model.ConstrainedSketch(
                        name='sketch_Pore_unvollkommen',
                        sheetSize=200,
                        transform=self.transform)
                    self.sketch_Pore_unvollkommen.Line(
                        point1=(-self.porenparameter_x/2.0*2.0, 0.0),
                        point2=(-self.porenparameter_x/2.0*3.0, 0.0))
                    self.sketch_Pore_unvollkommen.Line(
                        point1=(-self.porenparameter_x/2.0*3.0, 0.0),
                        point2=(-self.porenparameter_x/2.0*3.0, self.porenparameter_y/2.0*3.0))
                    self.sketch_Pore_unvollkommen.Line(
                        point1=(-self.porenparameter_x/2.0*3.0, self.porenparameter_y/2.0*3.0),
                        point2=(0.0, self.porenparameter_y/2.0*3.0))
                    self.sketch_Pore_unvollkommen.Line(
                        point1=(0.0, self.porenparameter_y/2.0*3.0),
                        point2=(0.0, self.porenparameter_y/2.0*2.0))
                    self.sketch_Pore_unvollkommen.Line(
                        point1=(0.0, self.porenparameter_y/2.0*2.0),
                        point2=(-self.porenparameter_x/2.0*2.0, self.porenparameter_y/2.0*2.0))
                    self.sketch_Pore_unvollkommen.Line(
                        point1=(-self.porenparameter_x/2.0*2.0, self.porenparameter_y/2.0*2.0),
                        point2=(-self.porenparameter_x/2.0*2.0, 0.0))
                    self.part_Pore_Skelett.CutExtrude(
                        sketchPlane=self.part_Pore_Skelett.faces[8],
                        sketchUpEdge=self.part_Pore_Skelett.edges[24],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_unvollkommen,
                        flipExtrudeDirection=OFF)
                    del self.sketch_Pore_unvollkommen
                    #Hilfsgeometrie loeschen (XZ-Ebene)
                    self.transform = self.part_Pore_Skelett.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett.faces[4],
                        sketchUpEdge=self.part_Pore_Skelett.edges[17],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_unvollkommen = model.ConstrainedSketch(
                        name='sketch_Pore_unvollkommen',
                        sheetSize=200,
                        transform=self.transform)
                    self.sketch_Pore_unvollkommen.rectangle(
                        point1=(0.0, self.porenparameter_z/2.0*2.0),
                        point2=(self.porenparameter_x/2.0*3.0, self.porenparameter_z/2.0*3.0))
                    self.part_Pore_Skelett.CutExtrude(
                        sketchPlane=self.part_Pore_Skelett.faces[4],
                        sketchUpEdge=self.part_Pore_Skelett.edges[17],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_unvollkommen,
                        flipExtrudeDirection=OFF)
                    del self.sketch_Pore_unvollkommen
                    #Ellipsoidskelett duplizieren
                    self.part_Pore_Skelett_2 = model.Part(
                        name='RVE_Pore_Skelett_2',
                        objectToCopy=self.part_Pore_Skelett,
                        compressFeatureList=ON)
                    #Solidifikation Skelett 1
                    self.part_Pore_Skelett.WirePolyLine(
                        points=((
                            self.part_Pore_Skelett.InterestingPoint(
                                edge=self.part_Pore_Skelett.edges[5],
                                rule=MIDDLE),
                            self.part_Pore_Skelett.InterestingPoint(
                                edge=self.part_Pore_Skelett.edges[6],
                                rule=MIDDLE)
                            ), ),
                        mergeType=IMPRINT,
                        meshable=ON)
                    self.part_Pore_Skelett.WirePolyLine(
                        points=((
                            self.part_Pore_Skelett.InterestingPoint(
                                edge=self.part_Pore_Skelett.edges[7],
                                rule=MIDDLE),
                            self.part_Pore_Skelett.InterestingPoint(
                                edge=self.part_Pore_Skelett.edges[10],
                                rule=MIDDLE)
                            ), ),
                        mergeType=IMPRINT,
                        meshable=ON)
                    self.part_Pore_Skelett.WirePolyLine(
                        points=((
                            self.part_Pore_Skelett.vertices[4],
                            self.part_Pore_Skelett.vertices[6]), ),
                        mergeType=IMPRINT,
                        meshable=ON)
                    self.part_Pore_Skelett.SolidLoft(
                        loftsections=(
                            (
                                self.part_Pore_Skelett.edges[8],
                                self.part_Pore_Skelett.edges[12],
                                self.part_Pore_Skelett.edges[14],
                                self.part_Pore_Skelett.edges[16],
                                self.part_Pore_Skelett.edges[17]),
                            (
                                self.part_Pore_Skelett.edges[0],
                                self.part_Pore_Skelett.edges[4],
                                self.part_Pore_Skelett.edges[15])),
                        paths=((self.part_Pore_Skelett.edges[3], ), ),
                        globalSmoothing=ON)
                    self.part_Pore_Skelett.PartitionCellByPlaneThreePoints(
                        point1=self.part_Pore_Skelett.vertices[0],
                        point3=self.part_Pore_Skelett.vertices[3],
                        cells=self.part_Pore_Skelett.cells,
                        point2=self.part_Pore_Skelett.InterestingPoint(
                            edge=self.part_Pore_Skelett.edges[11],
                            rule=MIDDLE))
                    self.transform = self.part_Pore_Skelett.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett.faces[0],
                        sketchUpEdge=self.part_Pore_Skelett.edges[11],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_Extrusion = model.ConstrainedSketch(
                        name='sketch_Pore_Extrusion',
                        sheetSize=200.0,
                        transform=self.transform)
                    del self.transform
                    self.sketch_Pore_Extrusion.rectangle(
                        point1=(0.0, 0.0),
                        point2=(-self.porenparameter_x/4.0, self.porenparameter_y/4.0))
                    self.part_Pore_Skelett.SolidExtrude(
                        sketchPlane=self.part_Pore_Skelett.faces[0],
                        sketchUpEdge=self.part_Pore_Skelett.edges[11],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_Extrusion,
                        depth=self.porenparameter_z/4.0,
                        flipExtrudeDirection=OFF)
                    del self.sketch_Pore_Extrusion
                    #Solidifikation Skelett 2
                    self.part_Pore_Skelett_2.WirePolyLine(
                        points=((
                            self.part_Pore_Skelett_2.InterestingPoint(
                                edge=self.part_Pore_Skelett_2.edges[4],
                                rule=MIDDLE),
                            self.part_Pore_Skelett.InterestingPoint(
                                edge=self.part_Pore_Skelett_2.edges[5],
                                rule=MIDDLE)
                            ), ),
                        mergeType=IMPRINT,
                        meshable=ON)
                    self.part_Pore_Skelett_2.WirePolyLine(
                        points=((
                            self.part_Pore_Skelett_2.InterestingPoint(
                                edge=self.part_Pore_Skelett_2.edges[9],
                                rule=MIDDLE),
                            self.part_Pore_Skelett_2.InterestingPoint(
                                edge=self.part_Pore_Skelett_2.edges[10],
                                rule=MIDDLE)
                            ), ),
                        mergeType=IMPRINT,
                        meshable=ON)
                    self.part_Pore_Skelett_2.WirePolyLine(
                        points=((
                            self.part_Pore_Skelett_2.vertices[2],
                            self.part_Pore_Skelett_2.vertices[5]), ),
                        mergeType=IMPRINT,
                        meshable=ON)
                    self.part_Pore_Skelett_2.SolidLoft(
                        loftsections=(
                            (
                                self.part_Pore_Skelett_2.edges[2],
                                self.part_Pore_Skelett_2.edges[6],
                                self.part_Pore_Skelett_2.edges[8],
                                self.part_Pore_Skelett_2.edges[14],
                                self.part_Pore_Skelett_2.edges[15]),
                            (
                                self.part_Pore_Skelett_2.edges[0],
                                self.part_Pore_Skelett_2.edges[12],
                                self.part_Pore_Skelett_2.edges[18])),
                        paths=((self.part_Pore_Skelett_2.edges[9], ), ),
                        globalSmoothing=ON)
                    self.part_Pore_Skelett_2.PartitionCellByPlaneThreePoints(
                        point1=self.part_Pore_Skelett_2.vertices[0],
                        point3=self.part_Pore_Skelett_2.vertices[3],
                        cells=self.part_Pore_Skelett_2.cells,
                        point2=self.part_Pore_Skelett_2.InterestingPoint(
                            edge=self.part_Pore_Skelett_2.edges[11],
                            rule=MIDDLE))
                    self.transform = self.part_Pore_Skelett_2.MakeSketchTransform(
                        sketchPlane=self.part_Pore_Skelett_2.faces[0],
                        sketchUpEdge=self.part_Pore_Skelett_2.edges[11],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_Extrusion = model.ConstrainedSketch(
                        name='sketch_Pore_Extrusion',
                        sheetSize=200.0,
                        transform=self.transform)
                    del self.transform
                    self.sketch_Pore_Extrusion.rectangle(
                        point1=(0.0, 0.0),
                        point2=(-self.porenparameter_z/4.0, self.porenparameter_y/4.0))
                    self.part_Pore_Skelett_2.SolidExtrude(
                        sketchPlane=self.part_Pore_Skelett_2.faces[0],
                        sketchUpEdge=self.part_Pore_Skelett_2.edges[11],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_Extrusion,
                        depth=self.porenparameter_x/4.0,
                        flipExtrudeDirection=ON)
                    del self.sketch_Pore_Extrusion
                    #2 Skelette der Ellipsoid fusionieren
                    self.assembly = model.rootAssembly
                    self.assembly.Instance(
                        name=self.name+'_Pore_Skelett_1',
                        part=self.part_Pore_Skelett,
                        dependent=ON)
                    self.assembly.Instance(
                        name=self.name+'_Pore_Skelett_2',
                        part=self.part_Pore_Skelett_2,
                        dependent=ON)
                    self.assembly.InstanceFromBooleanMerge(
                        name=self.name+'_Pore_unvollkommen',
                        instances=(
                            self.assembly.instances[self.name+'_Pore_Skelett_1'],
                            self.assembly.instances[self.name+'_Pore_Skelett_2'], ),
                        originalInstances=SUPPRESS,
                        domain=GEOMETRY)
                    self.part_Pore_unvollkommen = model.parts[self.name+'_Pore_unvollkommen']
                    self.part_Pore_unvollkommen.ReplaceFaces(
                        faceList = self.part_Pore_unvollkommen.faces[2:3]+self.part_Pore_unvollkommen.faces[4:5],
                        stitch=True)
                    self.assembly.deleteFeatures((
                        self.name+'_Pore_Skelett_1',
                        self.name+'_Pore_Skelett_2',
                        self.name+'_Pore_unvollkommen-1', ))
                    del model.parts[self.name+'_Pore_Skelett']
                    del model.parts[self.name+'_Pore_Skelett_2']
            elif (self.typ_Pore == 'Quader'):
                self.sketch_Pore.rectangle(
                    point1=(-self.porenparameter_x/2.0, -self.porenparameter_y/2.0),
                    point2=(self.porenparameter_x/2.0, self.porenparameter_y/2.0))
            elif (self.typ_Pore == 'Zylinder'):
                self.sketch_Pore.EllipseByCenterPerimeter(
                    center=(0.0, 0.0),
                    axisPoint1=(self.porenparameter_x/2.0, 0.0),
                    axisPoint2=(0.0, self.porenparameter_y/2.0))
            else:
                print('typ_Pore Error!')
            #Part Pore generieren
            self.part_Pore = model.Part(
                name=self.name+'_Pore',
                dimensionality=THREE_D,
                type=DEFORMABLE_BODY)
            if (self.typ_Pore == 'Ellipsoid' ):
                if (self.porenparameter_x == self.porenparameter_z):
                    self.part_Pore.BaseSolidRevolve(
                        sketch=self.sketch_Pore,
                        angle=360.0,
                        flipRevolveDirection=OFF)
                else:
                    del model.parts[self.name+'_Pore']
                    model.parts.changeKey(fromName=self.name+'_Pore_unvollkommen', toName=self.name+'_Pore')
                    self.part_Pore = model.parts[self.name+'_Pore']
                    #unvollkommene Teile loeschen
                    self.transform = self.part_Pore.MakeSketchTransform(
                        sketchPlane=self.part_Pore.faces[1],
                        sketchUpEdge=self.part_Pore.edges[19],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_cut = model.ConstrainedSketch(
                        name='sketch_Pore_Extrusion',
                        sheetSize=200.0,
                        transform=self.transform)
                    del self.transform
                    self.sketch_Pore_cut.rectangle(
                        point1=(0.0, 0.0),
                        point2=(-self.porenparameter_y/2.0, -self.porenparameter_x/2.0))
                    self.part_Pore.CutExtrude(
                        sketchPlane=self.part_Pore.faces[1],
                        sketchUpEdge=self.part_Pore.edges[19],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_cut,
                        flipExtrudeDirection=OFF)
                    del self.sketch_Pore_cut
                    self.transform = self.part_Pore.MakeSketchTransform(
                        sketchPlane=self.part_Pore.faces[0],
                        sketchUpEdge=self.part_Pore.edges[15],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        origin=(0.0, 0.0, 0.0))
                    self.sketch_Pore_cut = model.ConstrainedSketch(
                        name='sketch_Pore_Extrusion',
                        sheetSize=200.0,
                        transform=self.transform)
                    del self.transform
                    self.sketch_Pore_cut.rectangle(
                        point1=(0.0, self.porenparameter_z/2.0),
                        point2=(self.porenparameter_x/2.0, 0.0))
                    self.part_Pore.CutExtrude(
                        sketchPlane=self.part_Pore.faces[0],
                        sketchUpEdge=self.part_Pore.edges[15],
                        sketchPlaneSide=SIDE1,
                        sketchOrientation=RIGHT,
                        sketch=self.sketch_Pore_cut,
                        flipExtrudeDirection=OFF)
                    del self.sketch_Pore_cut
            elif (self.typ_Pore == 'Quader' or 'Zylinder'):
                self.part_Pore.BaseSolidExtrude(
                    sketch=self.sketch_Pore,
                    depth=self.porenparameter_z)
            #Assemble
            self.assembly = model.rootAssembly
            self.assembly.DatumCsysByDefault(CARTESIAN)
            self.assembly.Instance(
                name=self.name+'_Wuerfel',
                part=self.part_Wuerfel,
                dependent=ON)
            self.assembly.Instance(
                name=self.name+'_Pore',
                part=self.part_Pore,
                dependent=ON)
            if (self.typ_Pore == 'Ellipsoid' ):
                self.assembly.translate(
                    instanceList=(self.name+'_Pore', ),
                    vector=(0.0, 0.0, self.laenge_z/2.0))
            elif (self.typ_Pore == 'Quader' or 'Zylinder'):
                self.assembly.translate(
                    instanceList=(self.name+'_Pore', ),
                    vector=(0.0, 0.0, self.laenge_z/2.0-self.porenparameter_z/2.0))
            self.assembly.rotate(
                instanceList=(self.name+'_Pore', ),
                axisPoint=(0.0, 0.0, self.laenge_z/2.0),
                axisDirection=(1.0, 0.0, self.laenge_z/2.0),
                angle=self.porenparameter_rx)
            self.assembly.rotate(
                instanceList=(self.name+'_Pore', ),
                axisPoint=(0.0, 0.0, self.laenge_z/2.0),
                axisDirection=(0.0, 1.0, self.laenge_z/2.0),
                angle=self.porenparameter_ry)
            self.assembly.rotate(
                instanceList=(self.name+'_Pore', ),
                axisPoint=(0.0, 0.0, self.laenge_z/2.0),
                axisDirection=(0.0, 0.0, self.laenge_z/2.0+1),
                angle=self.porenparameter_rz)
            self.assembly.InstanceFromBooleanCut(
                name='RVE',
                instanceToBeCut=self.assembly.instances[self.name+'_Wuerfel'],
                cuttingInstances=(self.assembly.instances[self.name+'_Pore'], ),
                originalInstances=SUPPRESS)
            self.assembly.deleteFeatures((self.name+'_Wuerfel', self.name+'_Pore', ))
            del model.parts[self.name+'_Wuerfel']
            #del model.parts[self.name+'_Pore']
            self.assembly.translate(
                instanceList=(self.name+'-1', ),
                vector=(0.0, 0.0, -self.laenge_z/2.0))
            self.part_RVE = model.parts[self.name]
        elif (self.dimension == '2D'):
            #Sketch Wuerfel zeichnen
            self.sketch_Wuerfel = model.ConstrainedSketch(
                name='Seitenansicht_Wuerfel',
                sheetSize=200.0)
            self.sketch_Wuerfel.rectangle(
                point1=(0.0, 0.0),
                point2=(self.laenge_x/2.0, self.laenge_y/2.0)) #x- und y-Symmetrie
            #Part Wuerfel generieren
            self.part_Wuerfel = model.Part(
                name=self.name+'_Wuerfel',
                dimensionality=TWO_D_PLANAR,
                type=DEFORMABLE_BODY)
            self.part_Wuerfel.BaseShell(sketch=self.sketch_Wuerfel)
            #Sketch Pore zeichnen
            self.sketch_Pore = model.ConstrainedSketch(
                name='Seitenansicht_Pore',
                sheetSize=200.0)
            if (self.typ_Pore == 'Ellipsoid'):
                self.sketch_Pore.ConstructionLine(
                    point1=(0.0, -100.0),
                    point2=(0.0, 100.0))
                self.sketch_Pore.EllipseByCenterPerimeter(
                    center=(0.0, 0.0),
                    axisPoint1=(self.porenparameter_x/2.0, 0.0),
                    axisPoint2=(0.0, self.porenparameter_y/2.0))
                self.sketch_Pore.autoTrimCurve(
                    curve1=self.sketch_Pore.geometry[3],
                    point1=(-self.porenparameter_x/2.0, 0.0))
                self.sketch_Pore.Line(
                    point1=(0.0, self.porenparameter_y/2.0),
                    point2=(0.0, -self.porenparameter_y/2.0))
            elif (self.typ_Pore == 'Quader'):
                self.sketch_Pore.rectangle(
                    point1=(-self.porenparameter_x/2.0, -self.porenparameter_y/2.0),
                    point2=(self.porenparameter_x/2.0, self.porenparameter_y/2.0))
            elif (self.typ_Pore == 'Zylinder'):
                self.sketch_Pore.EllipseByCenterPerimeter(
                    center=(0.0, 0.0),
                    axisPoint1=(self.porenparameter_x/2.0, 0.0),
                    axisPoint2=(0.0, self.porenparameter_y/2.0))
            else:
                print('typ_Pore Error!')
            #Part Pore generieren
            self.part_Pore = model.Part(
                name=self.name+'_Pore',
                dimensionality=TWO_D_PLANAR,
                type=DEFORMABLE_BODY)
            self.part_Pore.BaseShell(sketch=self.sketch_Pore)
            #Assemble
            self.assembly = model.rootAssembly
            self.assembly.DatumCsysByDefault(CARTESIAN)
            self.assembly.Instance(
                name=self.name+'_Wuerfel',
                part=self.part_Wuerfel,
                dependent=ON)
            self.assembly.Instance(
                name=self.name+'_Pore',
                part=self.part_Pore,
                dependent=ON)
            self.assembly.rotate(
                instanceList=(self.name+'_Pore', ),
                axisPoint=(0.0, 0.0, self.laenge_z/2.0),
                axisDirection=(0.0, 0.0, self.laenge_z/2.0+1),
                angle=self.porenparameter_rz)
            self.assembly.InstanceFromBooleanCut(
                name='RVE',
                instanceToBeCut=self.assembly.instances[self.name+'_Wuerfel'],
                cuttingInstances=(self.assembly.instances[self.name+'_Pore'], ),
                originalInstances=SUPPRESS)
            self.assembly.deleteFeatures((self.name+'_Wuerfel', self.name+'_Pore', ))
            del model.parts[self.name+'_Wuerfel']
            #del model.parts[self.name+'_Pore']
            self.part_RVE = model.parts[self.name]
        else:
            print('dimension Error!')
    def set_und_surface(self):
        if (self.dimension == '3D'):
            self.part_RVE.Set(
                cells=self.part_RVE.cells.getSequenceFromMask(mask=('[#1 ]', ), ),
                name='Set_RVE')
        elif (self.dimension == '2D'):
            self.part_RVE.Set(
                faces=self.part_RVE.faces.getSequenceFromMask(mask=('[#1 ]', ), ),
                name='Set_RVE')
        else:
            print('dimension Error!')
    def vernetzen(self,global_Mesh_Size,poren_Mesh_Size):
        self.global_Mesh_Size = global_Mesh_Size
        self.poren_Mesh_Size = poren_Mesh_Size
        self.part_RVE.seedPart(
            size=self.global_Mesh_Size,
            deviationFactor=0.1,
            minSizeFactor=0.1)
        if (self.dimension == '3D'):
            if(self.typ_Pore == 'Ellipsoid'):
                self.part_RVE.seedEdgeBySize(
                    edges=self.part_RVE.edges.getSequenceFromMask(mask=('[#1ff ]', ), ),
                    size=self.poren_Mesh_Size,
                    deviationFactor=0.1,
                    minSizeFactor=0.1,
                    constraint=FINER)
            elif(self.typ_Pore == 'Quader'):
                if (self.porenparameter_rx == 0.0 and self.porenparameter_ry == 0.0 and self.porenparameter_rz == 0.0):
                    self.part_RVE.seedEdgeBySize(
                        edges=self.part_RVE.edges.getSequenceFromMask(mask=('[#1ff ]', ), ),
                        size=self.poren_Mesh_Size,
                        deviationFactor=0.1,
                        minSizeFactor=0.1,
                        constraint=FINER)
                else:
                    self.part_RVE.seedEdgeBySize(
                        edges=self.part_RVE.edges.getSequenceFromMask(mask=('[#3f ]', ), ),
                        size=self.poren_Mesh_Size,
                        deviationFactor=0.1,
                        minSizeFactor=0.1,
                        constraint=FINER)
            elif(self.typ_Pore == 'Zylinder'):
                self.part_RVE.seedEdgeBySize(
                    edges=self.part_RVE.edges.getSequenceFromMask(mask=('[#3f ]', ), ),
                    size=self.poren_Mesh_Size,
                    deviationFactor=0.1,
                    minSizeFactor=0.1,
                    constraint=FINER)
            self.part_RVE.setMeshControls(
                regions=self.part_RVE.cells,
                elemShape=TET,
                technique=FREE)
            self.elemType1 = ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
            self.elemType2 = ElemType(elemCode=C3D15, elemLibrary=STANDARD)
            self.elemType3 = ElemType(elemCode=C3D10, elemLibrary=STANDARD)
            self.part_RVE.setElementType(
                regions=self.part_RVE.sets['Set_RVE'],
                elemTypes=(self.elemType1, self.elemType2, self.elemType3))
            self.part_RVE.generateMesh()
        elif (self.dimension == '2D'):
            if(self.typ_Pore == 'Ellipsoid' or 'Zylinder'):
                self.part_RVE.seedEdgeBySize(
                    edges=self.part_RVE.edges.getSequenceFromMask(mask=('[#1 ]', ), ),
                    size=self.poren_Mesh_Size,
                    deviationFactor=0.1,
                    minSizeFactor=0.1,
                    constraint=FINER)
            elif(self.typ_Pore == 'Quader'):
                self.part_RVE.seedEdgeBySize(
                    edges=self.part_RVE.edges.getSequenceFromMask(mask=('[#3 ]', ), ),
                    size=self.poren_Mesh_Size,
                    deviationFactor=0.1,
                    minSizeFactor=0.1,
                    constraint=FINER)
            self.part_RVE.setMeshControls(
                regions=self.part_RVE.faces,
                    elemShape=TRI)
            self.elemType1 = ElemType(elemCode=CPE4R, elemLibrary=STANDARD)
            self.elemType2 = ElemType(elemCode=CPE3, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)
            self.part_RVE.setElementType(
                regions=self.part_RVE.sets['Set_RVE'],
                elemTypes=(self.elemType1, self.elemType2))
            self.part_RVE.generateMesh()

#Model defenieren
modelname = 'RVE'
RVE_global_Mesh_Size = 0.002
RVE_poren_Mesh_Size = 0.0004
ModelZuruecksetzen = 1

#-------------------------------------------------------------------------------
#Model zuruecksetzen
#-------------------------------------------------------------------------------
if (ModelZuruecksetzen == 0):
    if (modelname in mdb.models.keys()):
        del mdb.models[modelname]
    mdb.Model(name=modelname, modelType=STANDARD_EXPLICIT)
elif (ModelZuruecksetzen == 1):
    if (len(mdb.models.keys()) > 1):
        for i in range(0,len(mdb.models.keys())-1):
            del mdb.models[mdb.models.keys()[1]]
    mdb.models.changeKey(fromName=mdb.models.keys()[0], toName=modelname)

model = mdb.models[modelname]


#-------------------------------------------------------------------------------
#RVE-Geometrie
#-------------------------------------------------------------------------------
rve = RVE(
    name = 'RVE',
    dimension = '3D',
    laenge_x = 0.1,
    laenge_y = 0.1,
    laenge_z = 0.1,
    typ_Pore = 'Ellipsoid',
    porenparameter_x = 0.03,
    porenparameter_y = 0.03,
    porenparameter_z = 0.03,
    porenparameter_rx = 0.0,
    porenparameter_ry = 0.0,
    porenparameter_rz = 0.0)


rve.sketch_und_part()
rve.set_und_surface()
rve.vernetzen(
    global_Mesh_Size = RVE_global_Mesh_Size,
    poren_Mesh_Size = RVE_poren_Mesh_Size)
