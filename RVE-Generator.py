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

#Arbeitspfad definieren
arbeitspfad = r'W:\Hiwi\Li_Xinyang\GitHub\rve_generator'

#Arbeitspfad und Packages-Searching-Pfad aendern
os.chdir(arbeitspfad)
sys.path.insert(0, arbeitspfad+'/Packages')

#MATLAB vorbereiten
from mlab.releases import latest_release
from matlab import matlabroot
from mlab.releases import latest_release as matlab
matlab.addpath(matlab.genpath('MATLAB-Codes'))
matlab.addpath(matlab.genpath('Packages'))



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
    def sketch_und_part(self):
        if (self.dimension == '3D'):
            #Sketch Wuerfel zeichnen
            self.sketch_Wuerfel = model.ConstrainedSketch(
                name='Seitenansicht_Wuerfel',
                sheetSize=200.0)
            self.sketch_Wuerfel.rectangle(
                point1=(-self.laenge_x/2.0, -self.laenge_y/2.0),
                point2=(self.laenge_x/2.0, self.laenge_y/2.0))
            #Part Wuerfel generieren
            self.part_Wuerfel = model.Part(
                name=self.name+'_Wuerfel',
                dimensionality=THREE_D,
                type=DEFORMABLE_BODY)
            self.part_Wuerfel.BaseSolidExtrude(
                sketch=self.sketch_Wuerfel,
                depth=self.laenge_z/2.0) #z-Symmetrie
            #Sketch Pore zeichnen (fuer Quader und Zylinder)
            self.sketch_Pore = model.ConstrainedSketch(
                name='Seitenansicht_Pore',
                sheetSize=200.0)
            if (self.typ_Pore == 'Quader'):
                self.sketch_Pore.rectangle(
                    point1=(-self.porenparameter_x/2.0, -self.porenparameter_y/2.0),
                    point2=(self.porenparameter_x/2.0, self.porenparameter_y/2.0))
            elif (self.typ_Pore == 'Zylinder'):
                self.sketch_Pore.EllipseByCenterPerimeter(
                    center=(0.0, 0.0),
                    axisPoint1=(self.porenparameter_x/2.0, 0.0),
                    axisPoint2=(0.0, self.porenparameter_y/2.0))
            elif (self.typ_Pore == 'Ellipsoid' ):
                matlab.ellipsoidIgesOut(
                    self.porenparameter_x,
                    self.porenparameter_y,
                    self.porenparameter_z,
                    'Ellipsoid')
            #     if (self.porenparameter_x == self.porenparameter_z):
            #         self.sketch_Pore.ConstructionLine(
            #             point1=(0.0, -100.0),
            #             point2=(0.0, 100.0))
            #         self.sketch_Pore.EllipseByCenterPerimeter(
            #             center=(0.0, 0.0),
            #             axisPoint1=(self.porenparameter_x/2.0, 0.0),
            #             axisPoint2=(0.0, self.porenparameter_y/2.0))
            #         self.sketch_Pore.autoTrimCurve(
            #             curve1=self.sketch_Pore.geometry[3],
            #             point1=(-self.porenparameter_x/2.0, 0.0))
            #         self.sketch_Pore.Line(
            #             point1=(0.0, self.porenparameter_y/2.0),
            #             point2=(0.0, -self.porenparameter_y/2.0))
            else:
                print('typ_Pore Error!')
            #Part Pore generieren
            if (self.typ_Pore == 'Ellipsoid' ):
                # if (self.porenparameter_x == self.porenparameter_z):
                #     self.part_Pore.BaseSolidRevolve(
                #         sketch=self.sketch_Pore,
                #         angle=360.0,
                #         flipRevolveDirection=OFF)
                self.iges_Datei = mdb.openIges(
                    'Ellipsoid.igs',
                    msbo=False,
                    trimCurve=DEFAULT,
                    scaleFromFile=OFF)
                model.PartFromGeometryFile(
                    name=self.name+'_Pore',
                    geometryFile=self.iges_Datei,
                    combine=False,
                    stitchTolerance=1.0,
                    dimensionality=THREE_D,
                    type=DEFORMABLE_BODY,
                    convertToAnalytical=1,
                    stitchEdges=1,
                    scale=0.5) # Skalierung
                self.part_Pore = model.parts[self.name+'_Pore']
                self.part_Pore.AddCells(
                    faceList = self.part_Pore.faces,
                    flipped=False)
                del self.iges_Datei
                os.remove('abaqus_read_iges0.log') #Arbeitsordner aufraeumen
                os.remove('temp-Ellipsoid-new.sat')
                os.remove('Ellipsoid.igs')
            elif (self.typ_Pore == 'Quader' or 'Zylinder'):
                self.part_Pore = model.Part(
                    name=self.name+'_Pore',
                    dimensionality=THREE_D,
                    type=DEFORMABLE_BODY)
                self.part_Pore.BaseSolidExtrude(
                    sketch=self.sketch_Pore,
                    depth=self.porenparameter_z)
            #Assemble
            self.assembly = model.rootAssembly
            self.assembly.DatumCsysByDefault(CARTESIAN)
            self.assembly.Instance(
                name=self.name+'_Pore',
                part=self.part_Pore,
                dependent=ON)
            self.assembly.Instance(
                name=self.name+'_Wuerfel',
                part=self.part_Wuerfel,
                dependent=ON)
            #Translation
            self.assembly.translate(
                instanceList=(self.name+'_Wuerfel', ),
                vector=(0.0, 0.0, -self.laenge_z/2.0))
            if (self.typ_Pore == 'Ellipsoid'):
                self.assembly.translate(
                    instanceList=(self.name+'_Pore', ),
                    vector=(0.0, 0.0, 0.0))
            elif (self.typ_Pore == 'Quader' or 'Zylinder'):
                self.assembly.translate(
                    instanceList=(self.name+'_Pore', ),
                    vector=(0.0, 0.0, -self.porenparameter_z/2.0))
            #Rotation
            self.assembly.rotate(
                instanceList=(self.name+'_Pore', ),
                axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(1.0, 0.0, 0.0),
                angle=self.porenparameter_rx)
            self.assembly.rotate(
                instanceList=(self.name+'_Pore', ),
                axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(0.0, 1.0, 0.0),
                angle=self.porenparameter_ry)
            self.assembly.rotate(
                instanceList=(self.name+'_Pore', ),
                axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(0.0, 0.0,1.0),
                angle=self.porenparameter_rz)
            #Schneiden
            self.assembly.InstanceFromBooleanCut(
                name='RVE',
                instanceToBeCut=self.assembly.instances[self.name+'_Wuerfel'],
                cuttingInstances=(self.assembly.instances[self.name+'_Pore'], ),
                originalInstances=SUPPRESS)
            self.assembly.deleteFeatures((self.name+'_Wuerfel', self.name+'_Pore', ))
            # del model.parts[self.name+'_Wuerfel']
            # del model.parts[self.name+'_Pore']
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
                    edges=(
                        self.part_RVE.edges[0],
                        self.part_RVE.edges[1]),
                    size=self.poren_Mesh_Size,
                    deviationFactor=0.1,
                    minSizeFactor=0.1,
                    constraint=FINER)
            elif(self.typ_Pore == 'Quader'):
                self.part_RVE.seedEdgeBySize(
                    edges=self.part_RVE.edges.getSequenceFromMask(mask=('[#fff ]', ), ),
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
                technique=SWEEP)
            self.part_RVE.setSweepPath(region=self.part_RVE.cells[0], edge=self.part_RVE.edges[3], sense=REVERSE)
            self.elemType1 = ElemType(elemCode=C3D8T, elemLibrary=STANDARD)
            self.elemType2 = ElemType(elemCode=C3D6T, elemLibrary=STANDARD)
            self.elemType3 = ElemType(elemCode=C3D4T, elemLibrary=STANDARD)
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
RVE_global_Mesh_Size = 0.04
RVE_poren_Mesh_Size = 0.002



#-------------------------------------------------------------------------------
#RVE-Geometrie
#-------------------------------------------------------------------------------
rve = RVE(
    name = 'RVE',
    dimension = '3D',
    laenge_x = 0.4,
    laenge_y = 0.4,
    laenge_z = 0.4,
    typ_Pore = 'Ellipsoid',
    porenparameter_x = 0.06,
    porenparameter_y = 0.02,
    porenparameter_z = 0.04,
    porenparameter_rx = 0.0,
    porenparameter_ry = 0.0,
    porenparameter_rz = 45.0)

#-------------------------------------------------------------------------------
#Modelle zuruecksetzen
#-------------------------------------------------------------------------------
if (len(mdb.models.keys()) > 1):
    for i in range(0,len(mdb.models.keys())-1):
        del mdb.models[mdb.models.keys()[1]]

if (mdb.models.keys()[0] == modelname):
    mdb.models.changeKey(fromName=modelname, toName='Model-x')
    mdb.Model(name=modelname, modelType=STANDARD_EXPLICIT)
    del mdb.models['Model-x']
else:
    mdb.Model(name=modelname, modelType=STANDARD_EXPLICIT)
    del mdb.models[mdb.models.keys()[0]]

model = mdb.models[modelname]

#-------------------------------------------------------------------------------
#Preprocessing
#-------------------------------------------------------------------------------

rve.sketch_und_part()
rve.set_und_surface()
rve.vernetzen(
    global_Mesh_Size = RVE_global_Mesh_Size,
    poren_Mesh_Size = RVE_poren_Mesh_Size)
