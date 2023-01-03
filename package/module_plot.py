import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from .module_reso import *

def ploting():
    
    Valuex=[]
    Valuey=[]
    state=0
    def openNewWindow():
        
        def get_value(val):
            fig.clear()
            for j in range(0,int(np.size(Y[:,0])/3)) : 
                plt.axis('off')
                i=int(val)
                plt.plot(Y[j*3,:i], Y[j*3+1,:i],linewidth='0.75')
                plt.scatter(Y[j*3,i-1], Y[j*3+1,i-1])
                canvas.draw()
                
        fig = plt.figure()     
        plt.axis('off')
        for j in range(0,int(np.size(Y[:,0])/3)): 
            plt.plot(Y[j*3,:], Y[j*3+1,:],linewidth='0.75')
        newWindow = tkinter.Toplevel(root)
        newWindow.title("New Window")
        newWindow.resizable(width=False,height=False)

        F9 = tkinter.LabelFrame(newWindow, text = 'OBJ.1',labelanchor="n",padx =20)
        S1 = tkinter.Scale(newWindow, from_=0, to=np.size(Y[0,:]), orient="horizontal",length=300,command=get_value)
        canvas = FigureCanvasTkAgg(fig, master=F9)
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
                C1.create_rectangle(x-2,y-2,x+2,y+2,fill="red")
                if (showvalue==0):
                    C1.create_text(x,y-8,text="{},{}".format(x,400-y))
            else:
                C1.create_rectangle(x-2,y-2,x+2,y+2,fill="blue")
                C1.create_line(x,y,Valuex[-1],-(Valuey[-1]-400),fill="blue",width=1)
                if (showvalue==0):
                    C1.create_text(x,y-8,text="{},{}".format(x,400-y))
            Valuex.append(x)
            Valuey.append(400-y)
        def valider():
            Vit=np.array([[0,0,0]])
            Pos=np.array([[0,0,0]])
            statecheck1=varRK2.get()
            statecheck2=varRK4.get()
            statecheck3=VarAdaprkf45.get()
            stateMasses=VarMasses.get()
            nb_corps=int(len(Valuex)/2)
            for i in range(nb_corps): 
                    Vx=-Valuex[-(i*2+1)-1]+Valuex[-(i*2)-1]
                    Vy=-Valuey[-(i*2+1)-1]+Valuey[-(i*2)-1]
                    Vit=np.vstack((Vit,np.array([[Vx,Vy,0]])))
                    Pos=np.vstack((Pos,np.array([[Valuex[-(i*2+1)-1],Valuey[-(i*2+1)-1],0]])))
            Vit=Vit[1:,:]  
            Pos=Pos[1:,:]    
            y=np.vstack((Pos,Vit)).reshape(-1,1)
            N=int(E1.get()) 
            t_max=int(E2.get())
            Masses=[]
            if stateMasses==1:
                for i in range(nb_corps):
                    Masses.append(1000)
            else:        
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
            if (statecheck1==1 and statecheck2==statecheck3==0):
                Y=Rk2(y,t0=0,tf=t_max,N=N,nb_corps=nb_corps,M=Masses)  
            if (statecheck2==1 and statecheck1==statecheck3==0):
                Y=Rk4(y,t0=0,tf=t_max,N=N,nb_corps=nb_corps,M=Masses)    
            if (statecheck3==1 and statecheck1==statecheck2==0):
                Y=adaptive_rkf45(y,t0=0,tf=t_max,N=N,nb_corps=nb_corps,M=Masses)   
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
        E1.insert(0 , '5000')
        E2 = tkinter.Entry(F6)
        E2.insert(0 , '5')
        E4 = tkinter.Entry(F6)
        E4.insert(0 , '1000,1000,1000,1000')
        L25 = tkinter.Label(F6, text ='itération')
        L26 = tkinter.Label(F6, text = 'tf')
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
        L28.grid(row=2, column=0)
        E1.grid(row=0, column=1)
        E2.grid(row=1, column=1)
        E4.grid(row=2, column=1)
        F6.grid(row=0, column=1,columnspan=5)
        
        C1.pack()
        B1.grid(row=0, column=0)
        B2.grid(row=1, column=0)
        B3.grid(row=2, column=0)
        F4.grid(row=0, column=0)        
        F3.grid(row=2, column=0,columnspan=6)        
        varRK4 = tkinter.IntVar()
        varRK2 = tkinter.IntVar()
        VarAdaprkf45 = tkinter.IntVar()
        VarMasses = tkinter.IntVar()
        CB1= tkinter.Checkbutton(F6,text="rk2", variable=varRK2)
        CB1.grid(row=3, column=0)
        CB2= tkinter.Checkbutton(F6,text="Rk4", variable=varRK4)
        CB2.grid(row=3, column=1)
        CB3= tkinter.Checkbutton(F6,text="Rk45(adapt)", variable=VarAdaprkf45)
        CB3.grid(row=3, column=2)
        CB4= tkinter.Checkbutton(F6,text="Masses auto", variable=VarMasses)
        CB4.grid(row=2, column=2)
        
        def showvalue():
            for i in range(len(Valuex)):
                C1.create_text(Valuex[i],-(Valuey[i]-400)-8,text="{},{}".format(Valuex[i],-(Valuey[i]-400)))
        
        
        B4 = tkinter.Button(F4,text='Show value',command=showvalue)
        B4.grid(row=3, column=0)
                    
       
    def openNewWindowError():
        openNewWindowError = tkinter.Toplevel(root)
        openNewWindowError.title("Error")
        L1 = tkinter.Label(openNewWindowError, text = "Please check Rk6 ! OR ! Rk4")
        L1.pack()
        
    def openNewWindow3():
        def valider():  
            statecheck1=varRK4.get()
            statecheck2=varRK2.get()
            statecheck3=varRK45adap.get()
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
            if (statecheck1==1 and statecheck2==statecheck3==0):
                Y=Rk4(y,t0=0,tf=tf,N=N,nb_corps=nb_corps,M=M)  
            if (statecheck2==1 and statecheck1==statecheck3==0):
                Y=Rk2(y,t0=0,tf=tf,N=N,nb_corps=nb_corps,M=M)   
            if (statecheck3==1 and statecheck1==statecheck2==0):
                Y=adaptive_rkf45(y,t0=0,tf=tf,N=N,nb_corps=nb_corps,M=M)   
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
        F1.grid(row=1, column=0,columnspan=3)
        B1.grid(row=0, column=0)
        B2.grid(row=0, column=1)
        varRK4 = tkinter.IntVar()
        varRK2 = tkinter.IntVar()
        varRK45adap = tkinter.IntVar()
        CB1= tkinter.Checkbutton(newWindow3,text="rk2", variable=varRK2)
        CB1.grid(row=3, column=0)
        CB2= tkinter.Checkbutton(newWindow3,text="Rk4", variable=varRK4)
        CB2.grid(row=3, column=1)
        CB3= tkinter.Checkbutton(newWindow3,text="adaptive_rkf45", variable=varRK45adap)
        CB3.grid(row=3, column=2)
    root = tkinter.Tk()
    def quit_me():
            root.quit()
            root.destroy()
    root.protocol("WM_DELETE_WINDOW", quit_me)
    root.wm_title("Embedding in Tk")
    root.resizable(width=False,height=False)
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