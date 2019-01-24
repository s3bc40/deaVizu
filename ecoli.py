# Powered by Python 2.7

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts : 
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import tlp

# The updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# The pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# The runGraphScript(scriptFile, graph) function can be called to launch
# another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# The main(graph) function must be defined 
# to run the script on the current graph

#===================================
# PART I
#===================================
def preProcess(rootGraph,locus,label,size,color,regPositive,regNegative,layout,labelColor,labelBorderColor):
  """ Function for preprocessing the root graph before applying main algorithms.

  This function permits to get the first vizualisation of the root graph. It will apply
  labels corresponding to the locus at each nodes. A size to each nodes in order to be able to 
  see the labels. A color to each edges of the root graph in function of their regulation type.
  And finally, a layout algorithm, here the "Random Layout", to store a position for each nodes.

  Args:
    rootGraph (tlp.Graph) : the root graph of this vizualisation project
    locus (tlp.StringProperty) : a property linked to the loci of a gene, or a node
    label (tlp.StringProperty) : a property linked to the "viewLayout" of the root graph 
    size (tlp.DoubleProperty) : a property linked to the "viewSize" of the root graph 
    color (tlp.ColorProperty) : a property linked to the "viewColor" of the root graph 
    regPositive (tlp.BooleanProperty) : a property linked to the positive regulation of a gene
    regPositive (tlp.BooleanProperty) : a property linked to the negative regulation of a gene
    layout (tlp.LayoutProperty) : a property linked to the "viewLayout" of the root graph
    labelColor (tlp.ColorProperty) : a property linked to the "viewLabelColor" of the root graph 
    labelBorderColor (tlp.ColorProperty) : a property linked to the "viewLabelBorderColor" of the root graph

  Returns:
    None
  """
  for edge in rootGraph.getEdges():
    if(regNegative[edge] == True and regPositive[edge] == False):
      color[edge] = tlp.Color.Red
    elif(regNegative[edge] == False and regPositive[edge] == True):
      color[edge] = tlp.Color.Green
    elif(regNegative[edge] == True and regPositive[edge] == True):
      color[edge] = tlp.Color.Blue
    else:
      color[edge] = tlp.Color.Amber
  label.copy(locus)
  color.setAllNodeValue(tlp.Color.Gray)
  labelColor.setAllNodeValue(tlp.Color.White)
  labelBorderColor.setAllNodeValue(tlp.Color.White)
  size.setAllNodeValue(tlp.Size(5,5,0))
  rootGraph.applyLayoutAlgorithm("Random layout",layout)

#===================================
# PART II
#===================================

#===================================
# CONTRUCT TREE
def constructHierachicalTree(treeGraph,interactGraph):
  """ Function to initialize the recursive algorithm to construct the hierarchical tree with a root node.

  Args:
    treeGraph (tlp.Graph) : the hierarchical tree graph 
    interactGraph (tlp.Graph) : the gene interaction graph

  Returns:
    None
  """
  root = treeGraph.addNode()
  constructRecursiveTree(treeGraph,interactGraph,root)
  
def constructRecursiveTree(treeGraph,interactGraph,srcNode):
  """ Recursive algorithm to construct the tree graph by running through the gene interaction graph and its sugraphes.

  The principle of this algorithm is to construct a hierarchical tree where the nodes represents a cluster
  of genes and the leaves the nodes that we need to classify. Hence, the compute process consists to recursively
  roaming in the gene interactions subgraphes by adding a node in the tree graph for each subgraph encountered. When a subgraph
  doesn't have any child subgraph, all the nodes and edges are added and linked with the last cluster node in the tree graph.
  This algorithm end when all sugraphes from the gene interaction graph are computed.

  Args:
    treeGraph (tlp.Graph) : the hierarchical tree graph to construct
    interactGraph (tlp.Graph) : the gene interaction graph (or subgraph) to proceed down
    srcNode (tlp.Node) : the source node from the tree graph (initialized at root node)

  Returns:
    None
  """

  if(interactGraph.numberOfSubGraphs() == 0):
    for node in interactGraph.getNodes():
      treeGraph.addNode(node)
      treeGraph.addEdge(srcNode,node)
  else:
    for subGraph in interactGraph.getSubGraphs():
        childNode = treeGraph.addNode()
        treeGraph.addEdge(srcNode,childNode)
        constructRecursiveTree(treeGraph,subGraph,childNode)


def applyRadialAlgo(treeGraph,layout):
  """ Apply the radial algorithm available in Tulip.

  Args:
    treeGraph (tlp.Graph) : the hierarchical tree graph 
    layout (tlp.LayoutProperty) : a property linked to the "viewLayout" of the root graph

  Returns:
    None
  """
  treeGraph.applyLayoutAlgorithm("Tree Radial")

def colorNodes(graph,doubleMetric):
  colorProp = graph.getColorProperty("viewColor")
  colorScale = tlp.ColorScale([])
  colors = [tlp.Color.Green ,tlp.Color.Black,tlp.Color.Red]
  colorScale.setColorScale(colors)
  param = tlp.getDefaultPluginParameters("Color Mapping", graph)
  param["input property"] = doubleMetric
  param["color scale"] = colorScale
  graph.applyColorAlgorithm("Color Mapping", colorProp, param)

#===================================
# CONSTRUCT BUNDLES
def computeShortPath(treeGraph,srcNode,endNode):
  depth = treeGraph.getLocalIntegerProperty("depth");
  fromSrcNode = []
  fromEndNode = []
  computeShortPathRec(treeGraph,srcNode,endNode,fromSrcNode,fromEndNode,depth)
  return fromSrcNode + fromEndNode

def computeShortPathRec(treeGraph,srcNode,endNode,fromSrcNode,fromEndNode,depth):
  if(srcNode == endNode):
    fromSrcNode.append(srcNode)
    fromEndNode.reverse()
    return None
  depthSource = depth[srcNode]
  depthEnd = depth[endNode]
  if(depthSource > depthEnd):
    srcNode = treeGraph.getInNode(srcNode,1)
    fromSrcNode.append(srcNode)
  elif(depthSource < depthEnd):
    endNode = treeGraph.getInNode(endNode,1)
    fromEndNode.append(endNode)
  else:
    srcNode = treeGraph.getInNode(srcNode,1)
    fromSrcNode.append(srcNode)
    endNode = treeGraph.getInNode(endNode,1)
    fromEndNode.append(endNode)
  computeShortPathRec(treeGraph,srcNode,endNode,fromSrcNode,fromEndNode,depth)

def constructBundles(treeGraph,interactGraph,layout):
  for edge in interactGraph.getEdges():
    u,v = interactGraph.ends(edge)
    pathNode = computeShortPath(treeGraph,u,v)
    pathCoord = []
    for node in pathNode:
      pathCoord.append(layout[node])
    layout.setEdgeValue(edge,pathCoord)
      
#===================================
# PART III
#===================================
def createHierarchy(smallMultGraph,interactGraph,graph):
  listNodes = []
  for index in range(1,18):
    tp_i = graph.getLocalDoubleProperty("tp{} s".format(index))
    subCopyGraph = smallMultGraph.addSubGraph("tp{}".format(index))
    tlp.copyToGraph(subCopyGraph,interactGraph)
    metric = subCopyGraph.getLocalDoubleProperty("viewMetric")
    tp_i = graph.getDoubleProperty("tp{} s".format(index))
    metric.copy(tp_i)
    colorSmallMultiples(subCopyGraph,metric)

def colorSmallMultiples(graph,doubleMetric):
  colorProp = graph.getLocalColorProperty("viewColor")
  colorScale = tlp.ColorScale([])
  colors = [tlp.Color.Red, tlp.Color.Black ,tlp.Color.Green]
  colorScale.setColorScale(colors)
  param = tlp.getDefaultPluginParameters("Color Mapping", graph)
  param["input property"] = doubleMetric
  param["color scale"] = colorScale
  graph.applyColorAlgorithm("Color Mapping", colorProp, param)

def constructGrid(graph,columns):
	layout = graph.getLayoutProperty("viewLayout")
	bBox = tlp.computeBoundingBox(graph)
	numberSub = graph.numberOfSubGraphs()
	idSubgraph = 0
	line = 0
	while  idSubgraph < numberSub:
		for column in range(columns):
			if(idSubgraph < numberSub):
				idSubgraph += 1
				subGraph = graph.getSubGraph("tp{}".format(idSubgraph))
				drawSmallMultiple(subGraph,line,column,bBox,layout)
			else:
				break
		line += 1

def drawSmallMultiple(graph, line, column,bBox,layout):
	width = bBox.width()
	height = bBox.height()
	for node in graph.getNodes():
		layout[node] = tlp.Coord((column * width) + layout[node].getX(),(line * - height) + layout[node].getY(),0)
	for edge in graph.getEdges():
		newControlPoints = []
		for controlPoint in layout[edge]:
			controlPoint = tlp.Coord((column * width) + controlPoint.getX(),(line * - height) + controlPoint.getY(),0)
			newControlPoints.append(controlPoint)
		layout.setEdgeValue(edge,newControlPoints)

def createSmallMultiples(smallMultGraph, interactGraph):
  createHierarchy(smallMultGraph,interactGraph,graph)
  constructGrid(smallMultGraph,7)
  

#===================================
# MAIN
#===================================
def main(graph): 
  Locus = graph.getStringProperty("Locus")
  Negative = graph.getBooleanProperty("Negative")
  Positive = graph.getBooleanProperty("Positive")
  similarity = graph.getDoubleProperty("similarity")
  tp1_s = graph.getDoubleProperty("tp1 s")
  tp10_s = graph.getDoubleProperty("tp10 s")
  tp11_s = graph.getDoubleProperty("tp11 s")
  tp12_s = graph.getDoubleProperty("tp12 s")
  tp13_s = graph.getDoubleProperty("tp13 s")
  tp14_s = graph.getDoubleProperty("tp14 s")
  tp15_s = graph.getDoubleProperty("tp15 s")
  tp16_s = graph.getDoubleProperty("tp16 s")
  tp17_s = graph.getDoubleProperty("tp17 s")
  tp2_s = graph.getDoubleProperty("tp2 s")
  tp3_s = graph.getDoubleProperty("tp3 s")
  tp4_s = graph.getDoubleProperty("tp4 s")
  tp5_s = graph.getDoubleProperty("tp5 s")
  tp6_s = graph.getDoubleProperty("tp6 s")
  tp7_s = graph.getDoubleProperty("tp7 s")
  tp8_s = graph.getDoubleProperty("tp8 s")
  tp9_s = graph.getDoubleProperty("tp9 s")
  viewBorderColor = graph.getColorProperty("viewBorderColor")
  viewBorderWidth = graph.getDoubleProperty("viewBorderWidth")
  viewColor = graph.getColorProperty("viewColor")
  viewFont = graph.getStringProperty("viewFont")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewIcon = graph.getStringProperty("viewIcon")
  viewLabel = graph.getStringProperty("viewLabel")
  viewLabelBorderColor = graph.getColorProperty("viewLabelBorderColor")
  viewLabelBorderWidth = graph.getDoubleProperty("viewLabelBorderWidth")
  viewLabelColor = graph.getColorProperty("viewLabelColor")
  viewLabelPosition = graph.getIntegerProperty("viewLabelPosition")
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewMetric = graph.getDoubleProperty("viewMetric")
  viewRotation = graph.getDoubleProperty("viewRotation")
  viewSelection = graph.getBooleanProperty("viewSelection")
  viewShape = graph.getIntegerProperty("viewShape")
  viewSize = graph.getSizeProperty("viewSize")
  viewSrcAnchorShape = graph.getIntegerProperty("viewSrcAnchorShape")
  viewSrcAnchorSize = graph.getSizeProperty("viewSrcAnchorSize")
  viewTexture = graph.getStringProperty("viewTexture")
  viewTgtAnchorShape = graph.getIntegerProperty("viewTgtAnchorShape")
  viewTgtAnchorSize = graph.getSizeProperty("viewTgtAnchorSize")
  
  # Preprocess the data
  rootGraph = graph
  preProcess(rootGraph,Locus,viewLabel,viewSize,viewColor,Positive,Negative,viewLayout,viewLabelColor,viewLabelBorderColor)
  
  # Create Hierarchical Tree
  interactGraph = graph.getSubGraph("Genes interactions")
  treeGraph = graph.addSubGraph("Hierarchical Tree")
  constructHierachicalTree(treeGraph,interactGraph)
  applyRadialAlgo(treeGraph,viewLayout)
  colorNodes(treeGraph,tp1_s)

  # Construct bundles
  depth = treeGraph.getLocalIntegerProperty("depth");
  tlp.dagLevel(treeGraph,depth)
  constructBundles(treeGraph,interactGraph,viewLayout)

  # Create Small Multiple Graph
  smallMultGraph = graph.addSubGraph("Small Multiples")
  createSmallMultiples(smallMultGraph, interactGraph)
