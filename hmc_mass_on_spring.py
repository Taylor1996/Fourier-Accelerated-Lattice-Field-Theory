from vpython import*
#display(width=600,height=600,center=vector(6,0,0),background=color.white)
infile = open("ho_vis_output.txt", "r")
s = infile.read()
tokens = s.split('[')[1:]

# we discard tokens[0]
for i,path in enumerate(tokens):
    new_path = re.sub('\]|\n', '', path)
    tokens[i] = new_path
#x_data = [path for path in tokens[1:]]
pathList = [[0]*len(tokens[0])]*len(tokens)
for i in range(len(tokens)):
    pathList[i] = [float(x_i) for x_i in tokens[i].split()]

    

wall=box(pos=vector(0,1,0),size=vector(0.2,3,2),color=color.red)


t=0
dt=1

eq=vector(9,0,0)
pivot=vector(0,0,0)

for path in pathList:
    
    Mass=box(pos=vector(9+path[t],0,0),velocity=vector(0,0,0),size=vector(1,1,1),mass=1.0,color=color.blue)
    spring=helix(pos=pivot,axis=Mass.pos-pivot,radius=0.4,constant=1,thickness=0.1,coils=20,color=color.green)
    t+=dt
    while (t<1000):
      rate(100)
      #acc=(eq-Mass.pos)*(spring.constant/Mass.mass)
      #Mass.velocity=Mass.velocity+acc*dt
      Mass.pos= path[t]
      spring.axis=Mass.pos-spring.pos
      t=t+dt