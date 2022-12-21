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
            #if np.size(Valuex
            C1.create_rectangle(x,y,x+4,y+4,fill="red")
            Valuex.append(x)
            Valuey.append(400-y)
        def valider():
            statecheck1=varRK4.get()
            statecheck2=varRK2.get()
            if (np.size(Valuex)>1):
                L1.config(text = 'x = {}'.format(Valuex[0]))
                L2.config(text = 'y = {}'.format(Valuey[0]))
                Vx1=Valuex[1]-Valuex[0]
                L4.config(text = 'Vx = {}'.format(Vx1))
                Vy1=Valuey[1]-Valuey[0]
                L5.config(text = 'Vy = {}'.format(Vy1))
                Vit=np.array([[Vx1,Vy1,0]])
                Pos=np.array([[Valuex[0],Valuey[0],0]])
            if (np.size(Valuex)>3):
                L7.config(text = 'x = {}'.format(Valuex[2]))
                L8.config(text = 'y = {}'.format(Valuey[2]))
                Vx2=Valuex[3]-Valuex[2]
                L10.config(text = 'Vx = {}'.format(Vx2))
                Vy2=Valuey[3]-Valuey[2]
                L11.config(text = 'Vy = {}'.format(Vy2))
                Vit=np.vstack((Vit,np.array([[Vx2,Vy2,0]])))
                Pos=np.vstack((Pos,np.array([[Valuex[2],Valuey[2],0]])))
            if (np.size(Valuex)>5):
                L13.config(text = 'x = {}'.format(Valuex[4]))
                L14.config(text = 'y = {}'.format(Valuey[4]))
                Vx3=Valuex[5]-Valuex[4]
                L16.config(text = 'Vx = {}'.format(Vx3))
                Vy3=Valuey[5]-Valuey[4]
                L17.config(text = 'Vy = {}'.format(Vy3))
                Vit=np.vstack((Vit,np.array([[Vx3,Vy3,0]])))
                Pos=np.vstack((Pos,np.array([[Valuex[4],Valuey[4],0]])))
            if (np.size(Valuex)>7):
                L19.config(text = 'x = {}'.format(Valuex[6]))
                L20.config(text = 'y = {}'.format(Valuey[6]))                
                Vx4=Valuex[7]-Valuex[6]
                L22.config(text = 'Vx = {}'.format(Vx4))
                Vy4=Valuey[7]-Valuey[6]
                L23.config(text = 'Vy = {}'.format(Vy4))
                Vit=np.vstack((Vit,np.array([[Vx4,Vy4,0]])))
                Pos=np.vstack((Pos,np.array([[Valuex[6],Valuey[6],0]])))
                
            y=np.vstack((Pos,Vit)).reshape(-1,1)
            print(Pos)
            print(Vit)
            print(y)
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
        #root.geometry("415x500")
        
        
        F1 = tkinter.LabelFrame(newWindow2, text = 'OBJ.1',labelanchor="n",padx =20)
        L1 = tkinter.Label(F1, text = 'x = {}'.format(0))
        L2 = tkinter.Label(F1, text = 'y = {}'.format(0))
        L3 = tkinter.Label(F1, text = 'z = {}'.format(0))
        L4 = tkinter.Label(F1, text = 'Vx = {}'.format(0))
        L5 = tkinter.Label(F1, text = 'Vy = {}'.format(0))
        L6 = tkinter.Label(F1, text = 'Vz = {}'.format(0))
        
        F2 = tkinter.LabelFrame(newWindow2, text = 'OBJ.2',labelanchor="n",padx =20)
        L7 = tkinter.Label(F2, text = 'x = {}'.format(0))
        L8 = tkinter.Label(F2, text = 'y = {}'.format(0))
        L9 = tkinter.Label(F2, text = 'z = {}'.format(0))
        L10 = tkinter.Label(F2, text = 'Vx = {}'.format(0))
        L11 = tkinter.Label(F2, text = 'Vy = {}'.format(0))
        L12 = tkinter.Label(F2, text = 'Vz = {}'.format(0))
        
        F4 = tkinter.LabelFrame(newWindow2, text = 'OBJ.3',labelanchor="n",padx =20)
        L13 = tkinter.Label(F4, text = 'x = {}'.format(0))
        L14 = tkinter.Label(F4, text = 'y = {}'.format(0))
        L15 = tkinter.Label(F4, text = 'z = {}'.format(0))
        L16 = tkinter.Label(F4, text = 'Vx = {}'.format(0))
        L17 = tkinter.Label(F4, text = 'Vy = {}'.format(0))
        L18 = tkinter.Label(F4, text = 'Vz = {}'.format(0))
        
        F5 = tkinter.LabelFrame(newWindow2, text = 'OBJ.4',labelanchor="n",padx =20)
        L19 = tkinter.Label(F5, text = 'x = {}'.format(0))
        L20 = tkinter.Label(F5, text = 'y = {}'.format(0))
        L21 = tkinter.Label(F5, text = 'z = {}'.format(0))
        L22 = tkinter.Label(F5, text = 'Vx = {}'.format(0))
        L23 = tkinter.Label(F5, text = 'Vy = {}'.format(0))
        L24 = tkinter.Label(F5, text = 'Vz = {}'.format(0))
        
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
        B1 = tkinter.Button(newWindow2,text='Caluler',command=valider)
        B2 = tkinter.Button(newWindow2,text='Afficher',command=openNewWindow)
        C1 = tkinter.Canvas(F3,width=500,height=400,background="white")
        C1.bind("<Button-1>",actionclicsouris)

        L1.pack()
        L2.pack()
        L3.pack()
        L4.pack()
        L5.pack()
        L6.pack()
        F1.grid(row=1, column=0)
        
        L7.pack()
        L8.pack()
        L9.pack()
        L10.pack()
        L11.pack()
        L12.pack()
        F2.grid(row=1, column=1)
        
        L13.pack()
        L14.pack()
        L15.pack()
        L16.pack()
        L17.pack()
        L18.pack()
        F4.grid(row=1, column=2)
        
        L19.pack()
        L20.pack()
        L21.pack()
        L22.pack()
        L23.pack()
        L24.pack()
        F5.grid(row=1, column=3)
        
        L25.grid(row=1, column=0)
        L26.grid(row=2, column=0)
        L27.grid(row=3, column=0)
        L28.grid(row=4, column=0)
        E1.grid(row=1, column=1)
        E2.grid(row=2, column=1)
        E3.grid(row=3, column=1)
        E4.grid(row=4, column=1)
        F6.grid(row=1, column=4,columnspan=2)
        
        C1.pack()
        B1.grid(row=0, column=1)
        B2.grid(row=0, column=2)

        F3.grid(row=2, column=0,columnspan=6)        
        varRK4 = tkinter.IntVar()
        varRK2 = tkinter.IntVar()
        CB1= tkinter.Checkbutton(F6,text="rk4", variable=varRK4)
        CB1.grid(row=5, column=0)
        CB2= tkinter.Checkbutton(F6,text="rk2", variable=varRK2)
        CB2.grid(row=5, column=1)
        
        
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
    #root.geometry("415x500")
    def pickedfc(event):
        global valgg
        valgg=L1.get(L1.curselection())
        if (valgg=="SandBox"):
            B5.configure(command=openNewWindow2)  
        if (valgg=="RawInput"):    
            B5.configure(command=openNewWindow3)

    F1 = tkinter.LabelFrame(root, text = 'choose',labelanchor="n",padx =10,pady =10)
    L2 = tkinter.Label(root,text="Create by HUGO LANCERY")
    L1 = tkinter.Listbox(F1,height=2)
    L3 = tkinter.Listbox(F1,height=2)
    L1.bind('<<ListboxSelect>>',pickedfc)
    B5 = tkinter.Button(root,text='Continue')  
    L1.insert(1,"SandBox")
    L1.insert(2,"RawInput")
    L3.insert(1,"Rk2")
    L3.insert(2,"Rk4")
    F1.grid(row=1, column=0)
    L2.grid(row=0, column=0)
    L1.grid(row=0, column=0)
    L3.grid(row=0, column=1)
    B5.grid(row=2, column=0,columnspan=2,padx =10,pady =10)
    tkinter.mainloop()

    return

