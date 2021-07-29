import subprocess as sp
import io
import os
import shlex

class calc:
   '''
   calc -> this class implements attributes and methods that allow the user
      to run an electronic structure calculation with the code PW in the Quantum 
      Espresso package.

   * Please check https://www.quantum-espresso.org/Doc/INPUT_PW.html for a more
     complete description of PW input.

   Attributes:
   -----------

   control -> dictionary corresponding to &CONTROL section.
   system -> dictionary corresponding to &SYSTEM section.
   electrons -> dictionary corresponding to &ELECTRONS section.
   ions -> dictionary corresponding to &IONS section.
   cell -> dictionary corresponding to &CELL section.
   atomic_species -> list of the atomic species.
   atomic_mass -> list of the atomic mass of each species.
   pseudopotential -> list of the pseudopotential files.
   atom_type -> list associating each atom to its atomic species.
   x,y,z -> lists containing the initial cartesian positions of the atoms.
   if_pos1,if_pos2,if_pos3 -> lists of factors (either 1 or 0) to be multiplied 
      by the corresponding force component.
   kx,ky,xz -> lists containing the K vectors in the reciprocal space.
   wk -> list of the weights of the K vectors.
   v1,v2,v3 -> cell vectors.
   nk1,nk2,nk3 -> parameters specifying the K point grid in the Monhkhorst-Pack
      scheme.
   sk1,sk2,sk3 -> grid offsets; must be either 0 or 1.
   atomic_positions_units -> self-explanatory; check PW documentation for the 
      allowed values.
   cell_parameters_units -> self-explanatory; check PW documentation for the 
      allowed values.
   k_points_type -> self-explanatory; check PW documentation for the allowed 
      values.
      
   The attributes are initialized in the class constructor:
   '''
   def __init__(self):
      self.control=dict()
      self.system=dict()
      self.electrons=dict()
      self.ions=dict()
      self.cell=dict()
      self.atomic_species=[]
      self.atomic_mass=[]
      self.pseudopotential=[]
      self.atom_type=[]
      self.x=[]
      self.y=[]
      self.z=[]
      self.kx=[]
      self.ky=[]
      self.kz=[]
      self.wk=[]
      self.if_pos1=[]
      self.if_pos2=[]
      self.if_pos3=[]
      self.v1=[]
      self.v2=[]
      self.v3=[]
      self.nk1=None
      self.nk2=None
      self.nk3=None
      self.sk1=None
      self.sk2=None
      self.sk3=None
      self.atomic_positions_units=""
      self.cell_parameters_units=""
      self.k_points_type=""

   '''
   Public Methods:
   ---------------

   read_input(inpfile) -> reads PW input file 'inpfile' and assigns values to 
       attributes.
   build_input(saveinp,inpfile) -> returns PW input as a string. If 'saveinp' 
       is set to 'True', the input script will be saved into the file 'inpfile'.
   run(command_line,input_string,saveout,outfile,savecoords,coordfile) -> runs 
       'command_line', which should contain a path to PW executable, passing PW 
       input in 'input_string', returning an 'out' object. If 'saveout' is set 
       to 'True', the PW output will be dumped into the file 'outfile'. If 
       'savecoords' is set to 'True', last atomic coordinates in the PW output 
       will be saved into file 'coordfile'.
   '''   
   def build_input(self,saveinp=False,inpfile='input'):
      # Initialize the number of atoms and number of species
      self.system["nat"]=str(len(self.atom_type))
      self.system["ntyp"]=str(len(self.atomic_species))

      # Build &CONTROL section
      input_string="&CONTROL\n"

      for key in self.control:
         input_string+=key+" = "+str(self.control[key])+",\n"

      input_string+="/\n"

      # Build &SYSTEM section
      input_string+="&SYSTEM\n"

      for key in self.system:
         input_string+=key+" = "+str(self.system[key])+",\n"

      input_string+="/\n"

      # Build &ELECTRONS section
      input_string+="&ELECTRONS\n"

      for key in self.electrons:
         input_string+=key+" = "+str(self.electrons[key])+",\n"

      input_string+="/\n"

      # Build &IONS section
      input_string+="&IONS\n"

      for key in self.ions:
         input_string+=key+" = "+str(self.ions[key])+",\n"

      input_string+="/\n"

      # Build &CELL section
      if self.control["calculation"]=="'vc-relax'" or self.control["calculation"]=="'vc-md'":
         input_string+="&CELL\n"

         for key in self.cell:
            input_string+=key+" = "+str(self.cell[key])+",\n"

         input_string+="/\n"

      # Build card ATOMIC_SPECIES   
      input_string+="ATOMIC_SPECIES\n"
      
      for i in range(len(self.atomic_species)):
         input_string+=self.atomic_species[i]+"   "+str(self.atomic_mass[i])+"   "+self.pseudopotential[i]+"\n"

      # Build card CELL_PARAMETERS
      if int(self.system["ibrav"])==0:
         input_string+="CELL_PARAMETERS "+self.cell_parameters_units+"\n"
         input_string+=str(self.v1[0])+"   "+str(self.v1[1])+"   "+str(self.v1[2])+"\n"
         input_string+=str(self.v2[0])+"   "+str(self.v2[1])+"   "+str(self.v2[2])+"\n"
         input_string+=str(self.v3[0])+"   "+str(self.v3[1])+"   "+str(self.v3[2])+"\n"

      # Build card ATOMIC_POSITIONS
      input_string+="ATOMIC_POSITIONS "+self.atomic_positions_units+"\n"

      for i in range(len(self.atom_type)):
         if len(self.if_pos1)==0:
            input_string+=self.atom_type[i]+"   "+str(self.x[i])+"   "+str(self.y[i])+"   "+str(self.z[i])+"\n"
         else:
            input_string+=self.atom_type[i]+"   "+str(self.x[i])+"   "+str(self.y[i])+"   "+str(self.z[i])+"   "+\
            str(self.if_pos1[i])+"   "+str(self.if_pos2[i])+"   "+str(self.if_pos3[i])+"\n"

      # Build card K_POINTS
      input_string+="K_POINTS "+self.k_points_type+"\n"
      
      if self.k_points_type!="gamma":
         if self.k_points_type=="automatic":
            input_string+=str(self.nk1)+"   "+str(self.nk2)+"   "+str(self.nk3)+"   "+str(self.sk1)+"   "+\
            str(self.sk2)+"   "+str(self.sk3)+"\n"
         else:
            input_string+=str(self.wk)+"\n"

            for i in range(len(self.wk)):
               input_string+=str(self.kx[i])+"   "+str(self.ky[i])+"   "+str(self.kz[i])+"   "+str(self.wk[i])+"\n"

      if saveinp:
         self.__write_input__(inpfile,input_string)

      print("Done!")

      return input_string

   def read_input(self,inpfile="input"):
      iscontrol=False
      issystem=False
      iselectrons=False
      isions=False
      iscell=False
      iscellparameters=False
      isatomicspecies=False
      isatomicpositions=False
      iskpoints=False

      try:
         f=open(inpfile,"r")
      except FileNotFoundError:
         print("File not found!")

         return

      lines=f.readlines()

      j=0

      for line in lines:
         l=line.split()

         # Read &CONTROL section
         if l[0].lower()=="&control":	
            iscontrol=True
         elif iscontrol:
            if l[0]=="/":
               iscontrol=False
            else:
               self.control[l[0]]=l[2].replace(",","")

         # Read &SYSTEM section
         elif l[0].lower()=="&system":
            issystem=True
         elif issystem:
            if l[0]=="/":
               issystem=False
            else:
               self.system[l[0]]=l[2].replace(",","")

         # Read &ELECTRONS section
         elif l[0].lower()=="&electrons":
            iselectrons=True
         elif iselectrons:
            if l[0]=="/":
               iselectrons=False
            else:
               self.electrons[l[0]]=l[2].replace(",","")

         # Read &IONS section
         elif l[0].lower()=="&ions":
            isions=True
         elif isions:
            if l[0]=="/":
               isions=False
            else:
               self.ions[l[0]]=l[2].replace(",","")

         # Read &CELL section
         elif l[0].lower()=="&cell":
            iscell=True
         elif iscell:
            if l[0]=="/":
               iscell=False
            else:
               self.cell[l[0]]=l[2].replace(",","")

         # Read CELL_PARAMETERS card
         elif l[0].lower()=="cell_parameters":
            iscellparameters=True

            if len(l)==2:
               self.cell_parameters_units=l[1]

         elif iscellparameters:
            if j==0:            
               self.v1.append(float(l[0]))
               self.v1.append(float(l[1]))
               self.v1.append(float(l[2]))

               j=1            
            elif j==1:
               self.v2.append(float(l[0]))
               self.v2.append(float(l[1]))
               self.v2.append(float(l[2]))

               j=2
            elif j==2:
               self.v3.append(float(l[0]))
               self.v3.append(float(l[1]))
               self.v3.append(float(l[2]))

               j=0
               iscellparameters=False

         # Read ATOMIC_SPECIES card
         elif l[0].lower()=="atomic_species":
            isatomicspecies=True
         elif isatomicspecies:
            self.atomic_species.append(l[0])
            self.atomic_mass.append(l[1])
            self.pseudopotential.append(l[2])

            j+=1

            if j==int(self.system["ntyp"]):
               j=0
               isatomicspecies=False

         # Read ATOMIC_POSITIONS card
         elif l[0].lower()=="atomic_positions":
            isatomicpositions=True

            if len(l)==2:
               self.atomic_positions_units=l[1]

         elif isatomicpositions:
            self.atom_type.append(l[0])
            self.x.append(float(l[1]))
            self.y.append(float(l[2]))        
            self.z.append(float(l[3]))

            if len(l)==7:
               self.if_pos1.append(int(l[4])) 
               self.if_pos2.append(int(l[5])) 
               self.if_pos3.append(int(l[6]))

            j+=1

            if j==int(self.system["nat"]):
               j=0
               isatomicpositions=False

         # Read K_POINTS card
         elif l[0].lower()=="k_points":
            self.k_points_type=l[1]

            if self.k_points_type!="gamma":
               iskpoints=True
         elif iskpoints:
            if self.k_points_type=="automatic":
               self.nk1=int(l[0])
               self.nk2=int(l[1])
               self.nk3=int(l[2])
               self.sk1=int(l[3])
               self.sk2=int(l[4])
               self.sk3=int(l[5])
            else:
               if j==0:
                  nk=int(l[0])
                  j=1
               else:
                  self.kx.append(float(l[0]))
                  self.ky.append(float(l[1]))
                  self.kz.append(float(l[2]))
                  self.wk.append(float(l[3]))
   
                  j+=1

                  if j>nk:
                     iskpoints=False
                     j=0
   
      print("Done!")

      return

   def run(self,command_line,input_string,saveout=False,outfile="output",
           savecoords=False,coordfile="coord.xyz"):
      args=shlex.split(command_line)  
      pw=sp.Popen(args,stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE)

      try:
         input_bytes=bytes(input_string,"utf8")
      except:
          input_bytes=bytes(input_string)

      try:
         comm_pw=pw.communicate(input_bytes)
         out=io.StringIO(comm_pw[0].decode())
      except:
         raise PWException("Error when launching pw.x!")

      if os.path.isfile('CRASH'):
         raise PWException("pw.x has crashed!")
      else:
         print("Calculation finished!")

      pw_out_obj=self.__get_output_info__(out,saveout,outfile,savecoords,
                                          coordfile)

      return pw_out_obj

   '''
   Private methods:
   ----------------

   __write_input__(inpfile,input_string) -> writes PW input in 'input_string' 
       into file 'inpfile'.
   __get_output_info__(output,saveout,outfile,savecoords,coordfile) -> retrieves 
       information from the output of a PW calculation to be stored in an 'out' 
       object. If 'saveout' is set to 'True', the output will be saved into file 
       'outfile'. If 'savecoords' is set to 'True', last atomic coordinates will 
       be saved into 'coordfile'. 
   __write_coords__(coordfile,pw_out_obj) -> write last system coordinates taken from PW output stored 
      in the 'out' object 'pw_out_obj' into file 'outfile'.
   '''
   def __write_input__(self,inpfile,input_string):
      f=open(inpfile,"w")

      f.write(input_string)
      f.close()

      return

   def __get_output_info__(self,output,saveout=False,outfile="output",savecoords=False,coordfile="coord.xyz"):
      ry2ev=13.605
      bohr2angs=0.529177
      nextcellpar=False
      nextcoords=False
      nextforces=False
      i=0
      j=0
      a0=None
      t=[]
      x=[]
      y=[]
      z=[]
      fx=[]
      fy=[]
      fz=[]
      v1=[]
      v2=[]
      v3=[]
      natoms=None
      energy=[]
      magnetization=[]
      absmagnetization=[]
      efermi=None
      jobdone=False

      lines=output.readlines()

      if saveout:
         f=open(outfile,"w")

      for line in lines:
         if saveout:
            f.write(line)

         l=line.rstrip().split()

         if len(l)>=5 and l[0]=="!":
            energy.append(float(l[4])*ry2ev)
         elif len(l)>=5 and l[0].lower()=="the" and l[1].lower()=="fermi" and l[2].lower()=="energy":
            efermi=float(l[4])
         elif len(l)>=4 and l[1]=="magnetization":
            if l[0]=="total":
               try:
                  magnetization.append(float(l[3]))
               except:
                  pass
            elif l[0]=="absolute":
               try:
                  absmagnetization.append(float(l[3]))
               except:
                  pass
         elif len(l)>=5 and l[0].lower()=="lattice" and l[1].lower()=="parameter": 
            a0=float(l[4])*bohr2angs
         elif len(l)>=5 and l[0].lower()=="number" and l[1].lower()=="of":
            if l[2].lower()=="atoms/cell":
               natoms=int(l[4])           
         elif len(l)>=2 and l[0].lower()=="crystal" and l[1].lower()=="axes:":
            nextcellpar=True
            v1=[]
            v2=[]
            v3=[]
         elif nextcellpar:
            if i==0:
               v1.append(float(l[3])*a0)
               v1.append(float(l[4])*a0)
               v1.append(float(l[5])*a0)
            elif i==1:
               v2.append(float(l[3])*a0)
               v2.append(float(l[4])*a0)
               v2.append(float(l[5])*a0)
            elif i==2:
               v3.append(float(l[3])*a0)
               v3.append(float(l[4])*a0)
               v3.append(float(l[5])*a0)

               nextcellpar=False
            
            i+=1
         elif len(l)>=7 and l[0].lower()=="forces" and l[1].lower()=="acting" \
         and l[2].lower()=="on" and l[3].lower()=="atoms":
            nextforces=True
            fx=[]
            fy=[]
            fz=[]
            j=0
         elif nextforces:
            if len(l)>0:
               fx.append(float(l[6])*ry2ev/bohr2angs)
               fy.append(float(l[7])*ry2ev/bohr2angs)
               fz.append(float(l[8])*ry2ev/bohr2angs)

               j+=1

               if j==natoms:
                  j=0
                  nextforces=False
         elif len(l)>=1 and l[0].lower()=="atomic_positions":
            nextcoords=True
            t=[]
            x=[]
            y=[]
            z=[]
            j=0
         elif nextcoords:
            t.append(l[0])
            x.append(float(l[1]))
            y.append(float(l[2]))
            z.append(float(l[3]))
     
            j+=1
   
            if j==natoms:
               j=0
               nextcoords=False
         elif len(l)>=5 and l[3].lower()=="terminated":
             jobdone=True
   
      print("Data collected from PW output!")

      if saveout:
         f.close()

      pw_out_obj=out()
      pw_out_obj.atom_type=t
      pw_out_obj.x=x
      pw_out_obj.y=y
      pw_out_obj.z=z
      pw_out_obj.fx=fx
      pw_out_obj.fy=fy
      pw_out_obj.fz=fz
      pw_out_obj.v1=v1
      pw_out_obj.v2=v2
      pw_out_obj.v3=v3
      pw_out_obj.energy=energy
      pw_out_obj.magnetization=magnetization
      pw_out_obj.absmagnetization=absmagnetization
      pw_out_obj.efermi=efermi
      pw_out_obj.jobdone=jobdone

      if savecoords:
         self.__write_coords__(coordfile,pw_out_obj)

      return pw_out_obj

   def __write_coords__(self,coordfile,pw_out_obj):
      g=open(coordfile,"w")

      g.write("%d\n" % len(pw_out_obj.atom_type))
      g.write('Lattice="%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f" Properties=species:S:1:pos:R:3\n' \
      % (pw_out_obj.v1[0],pw_out_obj.v1[1],pw_out_obj.v1[2],pw_out_obj.v2[0],pw_out_obj.v2[1],\
      pw_out_obj.v2[2],pw_out_obj.v3[0],pw_out_obj.v3[1],pw_out_obj.v3[2]))

      for i in range(len(pw_out_obj.atom_type)):
         g.write("%s %.6f %.6f %.6f\n" % (pw_out_obj.atom_type[i],pw_out_obj.x[i],pw_out_obj.y[i],\
         pw_out_obj.z[i]))

      g.close()

      print("Coordinates saved!")

      return

class out:
   '''
   out -> This class stores in its attributes the values of variables taken
      from a PW calculation. An object of this class should be usually instantiated 
      from a call to 'calc.run()' method.

   Attributes:
   -----------

   atom_type -> list associating each atom to its atomic species.
   x,y,z -> lists containing the final cartesian positions of the atoms (Ang).
   fx,fy,fz -> list containing the final forces acting on atoms (ev/Ang).
   v1,v2,v3 -> cell vectors (units of lattice parameter).
   energy -> list containing the total energy of the system after every scf step (eV).
   efermi -> Fermi energy (eV).
   magnetization -> list containing the total magnetization of the system after 
      every scf step (Bohr magneton).
   absmagnetization -> list containing the absolute magnetization of the system
      after every scf step (Bohr magneton).
   jobdone -> a flag that is 'True' if the calculation finished properly or
      otherwise 'False'.   
      
   The attributes are initialized in the class constructor:   
   '''
   def __init__(self):
      self.atom_type=[]
      self.x=[]
      self.y=[]
      self.z=[]
      self.fx=[]
      self.fy=[]
      self.fz=[]
      self.v1=[]
      self.v2=[]
      self.v3=[]
      self.energy=[]
      self.efermi=None
      self.magnetization=[]
      self.absmagnetization=[]
      self.jobdone=False

class PWException(Exception):
   '''
   PWException -> defines an exception class for exeception handling.
   '''
   def __init__(self,message=""):
      self.message=message
      
   def __str__(self):
      return self.message
