class Mgrid():
    
    import numpy as np
    from scipy.sparse.linalg import LinearOperator
    
    def __init__(self,Nlevels,level = 1):
        self.Nlevels=Nlevels
        self.level = level
        self.c1 = 0.5625
        self.c2 = 0.1875
        self.c3 = 0.0625
        
    def apply_DBC(self,u):
        u[ 0,:] = -u[ 1,:]
        u[-1,:] = -u[-2,:]
        u[:, 0] = -u[:, 1]
        u[:,-1] = -u[:,-2]
        
    def jacobiRelaxation(self,nx,ny,u,f,iter=1,pre=False):
    
      dx=1.0/nx; dy=1.0/ny
      Ax=1.0/dx**2; Ay=1.0/dy**2
      Ap=1.0/(2.0*(Ax+Ay))

      #Dirichlet BC
      self.apply_DBC(u)

      if(pre and self.level>1):
        u[1:nx+1,1:ny+1] = -Ap*f[1:nx+1,1:ny+1]
        #Dirichlet BC
        self.apply_DBC(u)
        iter=iter-1

      for it in range(iter):
        u[1:nx+1,1:ny+1] = Ap*(Ax*(u[2:nx+2,1:ny+1] + u[0:nx,1:ny+1])
                             + Ay*(u[1:nx+1,2:ny+2] + u[1:nx+1,0:ny])
                             - f[1:nx+1,1:ny+1])
        #Dirichlet BC
        self.apply_DBC(u)

   

      res=np.zeros([nx+2,ny+2])
      res[1:nx+1,1:ny+1]=f[1:nx+1,1:ny+1]-(( Ax*(u[2:nx+2,1:ny+1]+u[0:nx,1:ny+1])
                                           + Ay*(u[1:nx+1,2:ny+2]+u[1:nx+1,0:ny])
                                           - 2.0*(Ax+Ay)*u[1:nx+1,1:ny+1]))
      return u,res
    
    def restriction(self,nx,ny,v):
        vc = np.zeros((nx+2,ny+2))
        
        vc[1:nx+1,1:ny+1]= 1/4*(v[1:2*nx:2,1:2*ny:2]\
                                         +v[1:2*nx:2,2:2*ny+1:2]+\
                                         v[2:2*nx+1:2,1:2*ny:2]+\
                                         v[2:2*nx+1:2,2:2*ny+1:2])
        
        return vc
      
    
    def prolongation(self,nx,ny,v):
        vt = np.zeros((2*nx+2,2*ny+2))
        
        vt[1:2*nx:2  ,1:2*ny:2] = self.c1*v[1:nx+1,1:ny+1]+\
                                            self.c2*(v[0:nx,1:ny+1]+v[1:nx+1,0:ny])\
                                            +self.c3*v[0:nx,0:ny]
        
        vt[2:2*nx+1:2,1:2*ny:2] = self.c1*v[1:nx+1,1:ny+1]+\
                                    self.c2*(v[2:nx+2,1:ny+1]+v[1:nx+1,0:ny])\
                                    +self.c3*v[2:nx+2,0:ny] 
    
    

        vt[1:2*nx:2  ,2:2*ny+1:2] = self.c1*v[1:nx+1,1:ny+1]+\
                                    self.c2*(v[0:nx  ,1:ny+1]+v[1:nx+1,2:ny+2])+\
                                    self.c3*v[0:nx  ,2:ny+2]
        
        vt[2:2*nx+1:2,2:2*ny+1:2] =self.c1*v[1:nx+1,1:ny+1]+\
                                    self.c2*(v[2:nx+2,1:ny+1]+v[1:nx+1,2:ny+2])+\
                                    self.c3*v[2:nx+2,2:ny+2]
        
        return vt
       
       
    
    def vycle(self,nx,ny,u,f):
        
        if (self.level == self.Nlevels):
            return self.jacobiRelaxation(nx,ny,u,f,iter = 50, pre = True)
   
        u, res = self.jacobiRelaxation(nx,ny,u,f,iter = 2, pre = True)
        
 

        res_coarse=self.restriction(nx//2,ny//2,res)

        error_coarse=np.zeros_like(res_coarse)
    
        self.level +=1

        error_coarse,res_coarse=self.vycle(nx//2,ny//2,error_coarse,res_coarse) 
    
        u+=self.prolongation(nx//2,ny//2,error_coarse)
         
        return  self.jacobiRelaxation(nx,ny,u,f,iter = 1, pre = False)


    