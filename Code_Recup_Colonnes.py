  dodo = graph.getDoubleProperty("tp17 s")
  liste = []
  for nodes in graph.getNodes():
    liste.append(dodo[nodes])
  with(open("test1.txt","w")) as file:
    for i in range (len(liste)):
      file.write(str(liste[i]) + "\n")
 # for i in range (len(liste)):
  #  print (liste[i])
