import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import customtkinter
import tkinter
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


def force(m1,X1,m2,X2) :
    G=4*np.pi**2    
    d=np.linalg.norm(X1-X2)
    f=-(G*m1*m2/(d**3))*(X1-X2)
    return f
    
def f(t,y,M) :
    nb_corps=len(M)
    F=np.zeros((nb_corps*6,1))
    for i in range(nb_corps) :
        S=np.zeros((3,1))
        for j in range(nb_corps) : 
            #print(y[i*3:(i+1)*3])
            if j != i:
                S+=force(M[i], y[i*3:(i+1)*3,0].reshape(-1,1),M[j],y[j*3:(j+1)*3,0].reshape(-1,1))
            
        F[i*3:(i+1)*3,0]=y[nb_corps*3+i*3:(i+1)*3+nb_corps*3,0]  
        F[nb_corps*3+i*3:(i+1)*3+nb_corps*3,0]=(1/M[i]*S).reshape(-1)    
    return F    

def Rk4(y,t0,tf,N,nb_corps,M):
    h=(tf-t0)/N
    t=t0
    Y=np.zeros((3*nb_corps,N))
    for i in range(0,N) : 
        k1=h*f(t,y,M)
        k2=h*f(t+h/2, y+k1/2,M)
        k3=h*f(t+h/2, y+k2/2,M)
        k4=h*f(t+h, y+k3,M)
        y=y+(1/6)*(k1 +2*k2 +2*k3 + k4)
        t+=h    
        for j in range(0,nb_corps) :
            Y[j*3:(j+1)*3,i]=y[j*3:(j+1)*3,0]
    return(Y)   

def Rk2(y,t0,tf,N,nb_corps,M):
    h=(tf-t0)/N
    t=t0
    Y=np.zeros((3*nb_corps,N))
    for i in range(0,N) : 
        k1=h*f(t,y,M)
        k2=h*f(t+h, y+k1,M)
        y=y+(k1+k2)/2
        t+=h    
        for j in range(0,nb_corps) :
            Y[j*3:(j+1)*3,i]=y[j*3:(j+1)*3,0]
    return(Y)  

def Verlet(y,t0,tf,N,nb_corps,M):
    h=(tf-t0)/N
    t=t0
    Y=np.zeros((3*nb_corps,N))
    for i in range(0,N) : 
        k1=h*f(t,y,M)
        k2=h*f(t+h, y+k1,M)
        y=y+(k1+k2)/2
        t+=h    
        for j in range(0,nb_corps) :
            Y[j*3:(j+1)*3,i]=y[j*3:(j+1)*3,0]
    return(Y)  


def ploting():
    
    Valuex=[]
    Valuey=[]
    state=0
    def openNewWindow():
        
        def get_value(val):
            fig.clear()
            for j in range(0,int(np.size(Y[:,0])/3)) : 
                i=int(val)
                ##plt.scatter(Y[j*3,i], Y[j*3+1,i], c=dico[j])
                plt.plot(Y[j*3,:i], Y[j*3+1,:i],c=dico[j],linewidth='0.75')
                canvas.draw()
                #print(val)
                
        fig = plt.figure()      
        #ax = plt.axes(projection='3d')
        dico={0:"y",1:"b",2:"r",3:"c",4:"m",5:"g"}
        #plt.axis('off')
        for j in range(0,int(np.size(Y[:,0])/3)): 
            plt.plot(Y[j*3,:], Y[j*3+1,:],c=dico[j],linewidth='0.75')
        
        newWindow = tkinter.Toplevel(root)
        newWindow.title("New Window")
        newWindow.resizable(width=False,height=False)

        F9 = tkinter.LabelFrame(newWindow, text = 'OBJ.1',labelanchor="n",padx =20)
        S1 = tkinter.Scale(newWindow, from_=0, to=np.size(Y[0,:]), orient="horizontal",length=300,command=get_value)
        canvas = FigureCanvasTkAgg(fig, master=F9)  # A tk.DrawingArea.
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, F9)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        F9.grid(row=0, column=0)
        S1.grid(row=1, column=0)
        
    def openNewWindow2():
        def actionclicsouris(event):
            
            C1.focus_set()
            x=event.x
            y=event.y
            if (np.size(Valuex)%2==0):
                C1.create_rectangle(x,y,x+4,y+4,fill="red")
                if (showvalue==0):
                    C1.create_text(x,y-8,text="{},{}".format(x,400-y))
            else:
                C1.create_rectangle(x,y,x+4,y+4,fill="blue")
                C1.create_line(x,y,Valuex[-1],-(Valuey[-1]-400),fill="blue",width=1)
                h=(int(E2.get())-0)/int(E1.get())
                v=np.round((((x-Valuex[-1])**2+(400-y-Valuey[-1])**2)**0.5)*h,2)
                if (showvalue==0):
                    C1.create_text(x,y-8,text="{},{}".format(x,400-y))
            Valuex.append(x)
            Valuey.append(400-y)
        def valider():
            Vit=np.array([[0,0,0]])
            Pos=np.array([[0,0,0]])
            statecheck1=varRK4.get()
            statecheck2=varRK2.get()
            nb_corps=int(E3.get())
            for i in range(nb_corps): 
                Vx=Valuex[-1-(i*2+1)]-Valuex[-1-(i*2)]
                Vy=Valuey[-1-(i*2+1)]-Valuey[-1-(i*2)]
                Vit=np.vstack((Vit,np.array([[Vx,Vy,0]])))
                Pos=np.vstack((Pos,np.array([[Valuex[-1-(i*2)],Valuey[-1-(i*2)],0]])))
            print(Vit)
            Vit=Vit[1:,:]        
            Pos=Pos[1:,:]    
            print(Vit) 
            y=np.vstack((Pos,Vit)).reshape(-1,1)
            N=int(E1.get())
            tf=int(E2.get())
            nb_corps=int(E3.get())
            Masses=[]
            mass=E4.get()
            mass=mass+","
            valinte=""            
            for i in mass:
                if (i!=","):
                    valinte=valinte+i
                if (i==","):
                    Masses.append(int(valinte))
                    valinte="" 
            global Y 
            if (statecheck1==1 and statecheck2==0):
                Y=Rk4(y,t0=0,tf=tf,N=N,nb_corps=nb_corps,M=Masses)  
            if (statecheck1==0 and statecheck2==1):
                Y=Rk2(y,t0=0,tf=tf,N=N,nb_corps=nb_corps,M=Masses)   
            print("done")
        newWindow2 = tkinter.Toplevel(root)
        newWindow2.title("New Window2")
        
        newWindow2.resizable(width=False,height=False)
        
        def reset():
            newWindow2.destroy()
            openNewWindow2()
            return
        
        F6 = tkinter.LabelFrame(newWindow2, text = 'Simulation',labelanchor="n",padx =20,pady =8.5)
        E1 = tkinter.Entry(F6)
        E1.insert(0 , '2000')
        E2 = tkinter.Entry(F6)
        E2.insert(0 , '3')
        E3 = tkinter.Entry(F6)
        E3.insert(0 , '4')
        E4 = tkinter.Entry(F6)
        E4.insert(0 , '10,10,10,10')
        L25 = tkinter.Label(F6, text ='itération')
        L26 = tkinter.Label(F6, text = 'tf')
        L27 = tkinter.Label(F6, text = 'nb Obj')
        L28 = tkinter.Label(F6, text = 'masses')
  
        F3 = tkinter.LabelFrame(newWindow2,padx =10,pady=10)
        F4 = tkinter.LabelFrame(newWindow2,padx =10,pady=10)
        B1 = tkinter.Button(F4,text='Caluler',command=valider)
        B2 = tkinter.Button(F4,text='Afficher',command=openNewWindow)
        B3 = tkinter.Button(F4,text='Reset',command=reset)
        C1 = tkinter.Canvas(F3,width=500,height=400,background="white")
        C1.bind("<Button-1>",actionclicsouris)

       
        
        L25.grid(row=0, column=0)
        L26.grid(row=1, column=0)
        L27.grid(row=2, column=0)
        L28.grid(row=3, column=0)
        E1.grid(row=0, column=1)
        E2.grid(row=1, column=1)
        E3.grid(row=2, column=1)
        E4.grid(row=3, column=1)
        F6.grid(row=0, column=1,columnspan=3)
        
        C1.pack()
        B1.grid(row=0, column=0)
        B2.grid(row=1, column=0)
        B3.grid(row=2, column=0)
        F4.grid(row=0, column=0)        
        F3.grid(row=2, column=0,columnspan=6)        
        varRK4 = tkinter.IntVar()
        varRK2 = tkinter.IntVar()
        CB1= tkinter.Checkbutton(F6,text="rk4", variable=varRK4)
        CB1.grid(row=5, column=0)
        CB2= tkinter.Checkbutton(F6,text="rk2", variable=varRK2)
        CB2.grid(row=5, column=1)
        
        def showvalue():
            for i in range(len(Valuex)):
                C1.create_text(Valuex[i],-(Valuey[i]-400)-8,text="{},{}".format(Valuex[i],-(Valuey[i]-400)))
        
        
        B4 = tkinter.Button(F4,text='Show value',command=showvalue)
        B4.grid(row=3, column=0)
        #h=(int(E2.get())-0)/int(E1.get())
        #v=np.round((((Valuex+(Valuex[-1]-400))**2+(400-y-Valuey[-1])**2)**0.5)*h,2)
                    
       
    def openNewWindowError():
        openNewWindowError = tkinter.Toplevel(root)
        openNewWindowError.title("Error")
        L1 = tkinter.Label(openNewWindowError, text = "Please check Rk2 ! OR ! Rk4")
        L1.pack()
        
    def openNewWindow3():
        def valider():  
            statecheck1=varRK4.get()
            statecheck2=varRK2.get()
            Pos=[]  
            Vit=[]
            M=[]
            PosStr=E1.get()
            PosStr=PosStr+","
            VitStr=E2.get()
            VitStr=VitStr+","
            MStr=E6.get()
            MStr=MStr+","
            N=int(E3.get())
            global nb_corps
            tf=int(E4.get())
            nb_corps=int(E5.get())
            
            valinte=""
            valinte1=""
            valinte2=""            
            for i in PosStr:
                if (i!=","):
                    valinte2=valinte2+i
                if (i==","):
                    Pos.append(int(valinte2))
                    valinte2="" 
            for j in VitStr:
                if (j!=","):
                    valinte1=valinte1+j
                if (j==","):
                    Vit.append(int(valinte1)) 
                    valinte1="" 
            for k in MStr:
                if (k!=","):
                    valinte=valinte+k
                if (k==","):
                    M.append(int(valinte)) 
                    valinte=""         
            y=np.vstack((Pos,Vit)).reshape(-1,1)
            global Y 
            if (statecheck1==1 and statecheck2==0):
                Y=Rk4(y,t0=0,tf=tf,N=N,nb_corps=nb_corps,M=M)  
            if (statecheck1==0 and statecheck2==1):
                Y=Rk2(y,t0=0,tf=tf,N=N,nb_corps=nb_corps,M=M)  
        newWindow3 = tkinter.Toplevel(root)
        newWindow3.title("New Window")
        F1 = tkinter.LabelFrame(newWindow3, text = 'OBJ.i',labelanchor="n",padx =10,pady =10)
        E1 = tkinter.Entry(F1)
        E1.insert(0 , '2,5,0,5,2,0,3,3,0,1,3,0')
        E2 = tkinter.Entry(F1)
        E2.insert(0 , '5,5,0,5,2,0,1,5,0,5,5,0')
        E3 = tkinter.Entry(F1)
        E3.insert(0 , '1000')
        E4 = tkinter.Entry(F1)
        E4.insert(0 , '2')
        E5 = tkinter.Entry(F1)
        E5.insert(0 , '4')
        E6 = tkinter.Entry(F1)
        E6.insert(0 , '1,1,1,1')
        L1 = tkinter.Label(F1, text = 'xi,yi,zi')
        L2 = tkinter.Label(F1, text = 'Vxi,Vyi,Vzi')
        L3 = tkinter.Label(F1, text ='itération')
        L4 = tkinter.Label(F1, text = 'tf')
        L5 = tkinter.Label(F1, text = 'nb Obj')
        L6 = tkinter.Label(F1, text = 'masses')
        L7 = tkinter.Label(F1, text = 'format : [x1,y1,z1,x2,...]') 
        B1 = tkinter.Button(newWindow3,text='Caluler',command=valider)
        B2 = tkinter.Button(newWindow3,text='Afficher',command=openNewWindow)

        L7.grid(row=0, column=0,columnspan=2)
        L1.grid(row=1, column=0)
        L2.grid(row=2, column=0)
        L3.grid(row=3, column=0)
        L4.grid(row=4, column=0)
        L5.grid(row=5, column=0)
        L6.grid(row=6, column=0)
        E1.grid(row=1, column=1)
        E2.grid(row=2, column=1)
        E3.grid(row=3, column=1)
        E4.grid(row=4, column=1)
        E5.grid(row=5, column=1)
        E6.grid(row=6, column=1)
        F1.grid(row=1, column=0,columnspan=2)
        B1.grid(row=0, column=0)
        B2.grid(row=0, column=1)
        varRK4 = tkinter.IntVar()
        varRK2 = tkinter.IntVar()
        CB1= tkinter.Checkbutton(newWindow3,text="rk4", variable=varRK4)
        CB1.grid(row=3, column=0)
        CB2= tkinter.Checkbutton(newWindow3,text="rk2", variable=varRK2)
        CB2.grid(row=3, column=1)
    root = tkinter.Tk()
    def quit_me():
            root.quit()
            root.destroy()
    root.protocol("WM_DELETE_WINDOW", quit_me)
    root.wm_title("Embedding in Tk")
    root.resizable(width=False,height=False)
    #root.geometry("415x500")
    def pickedfc(event):
        global valgg
        valgg=L1.get(L1.curselection())
        if (valgg=="SandBox"):
            B5.configure(command=openNewWindow2)  
        if (valgg=="RawInput"):    
            B5.configure(command=openNewWindow3)

    F1 = tkinter.LabelFrame(root, text = 'choose',labelanchor="n",padx =10,pady =10)
    L2 = tkinter.Label(root,text="N trajectory model",padx =10,pady =10)
    L1 = tkinter.Listbox(F1,height=2)
    L1.bind('<<ListboxSelect>>',pickedfc)
    B5 = tkinter.Button(root,text='Continue')  
    L1.insert(1,"SandBox")
    L1.insert(2,"RawInput")
    F1.grid(row=1, column=0)
    L2.grid(row=0, column=0)
    L1.grid(row=0, column=0)
    B5.grid(row=2, column=0,columnspan=2,padx =10,pady =10)
    tkinter.mainloop()

    return

ploting()
